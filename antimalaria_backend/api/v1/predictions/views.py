from django.http import HttpResponse
from rest_framework import status, viewsets
from rest_framework.views import APIView
from rest_framework.decorators import action
import csv
import io
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django.conf import settings
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Crippen
from rest_framework.response import Response
from concurrent.futures import ThreadPoolExecutor, as_completed
from .serializers import PredictionSerializer
from rest_framework.permissions import IsAuthenticated
from api.models import Prediction, Compound, PredictionCompound, MLModel
from drf_spectacular.utils import extend_schema_view, extend_schema
from .schemas import (
    prediction_list_schema,
    prediction_retrieve_schema,
    prediction_destroy_schema,
    predict_ic50_schema,
    prediction_download_schema
)

@extend_schema_view(
    list=extend_schema(**prediction_list_schema),
    retrieve=extend_schema(**prediction_retrieve_schema),
    destroy=extend_schema(**prediction_destroy_schema)
)
class PredictionViewSet(viewsets.ModelViewSet):
    serializer_class = PredictionSerializer
    permission_classes = [IsAuthenticated]
    http_method_names = ['get', 'delete', 'head', 'options']

    def get_queryset(self):
        qs = Prediction.objects.all() if self.request.user.role == 'admin' else Prediction.objects.filter(user=self.request.user)
        return qs.select_related('ml_model', 'user').prefetch_related('prediction_compounds__compound').order_by('-created_at') 

    def list(self, request, *args, **kwargs):
        queryset = self.filter_queryset(self.get_queryset())
        serializer = self.get_serializer(queryset, many=True)
        return Response({
            "status": "success",
            "message": "Predictions retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)
    
    def retrieve(self, request, *args, **kwargs):
        instance = self.get_object()
        serializer = self.get_serializer(instance)
        return Response({
            "status": "success",
            "message": "Prediction retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)
    
    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        instance.delete()
        return Response({
            "status": "success",
            "message": "Prediction deleted successfully."
        }, status=status.HTTP_204_NO_CONTENT)
    
    
    @action(detail=False, methods=["delete"], url_path="delete-all")
    def delete_all(self, request):
        user = request.user
        qs = Prediction.objects.filter(user=user)
        count, _ = qs.delete()

        return Response(
            {
                "status": "success",
                "message": f"All your predictions have been deleted ({count} entries)."
            },
            status=200
        )

class PredictIC50View(APIView):
    permission_classes = [IsAuthenticated]

    @extend_schema(
        summary="Predict IC50 values for SMILES",
        description="Predict IC50 for valid SMILES and mark invalid ones as NaN, but still store them in DB.",
    )
    def post(self, request, *args, **kwargs):
        user = request.user
        csv_file = request.FILES.get("file")
        smiles_input = request.data.get("smiles")
        model_descriptor = request.data.get("model_descriptor")
        model_method = request.data.get("model_method")

        if not model_descriptor or not model_method:
            return Response(
                {"error": "model_descriptor and model_method are required."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # 🔹 Ambil SMILES dari CSV atau teks
        smiles_list = []
        if csv_file:
            if not csv_file.name.endswith(".csv"):
                return Response(
                    {"error": "Only CSV files are supported."},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            try:
                decoded_file = csv_file.read().decode("utf-8-sig")
                reader = csv.reader(io.StringIO(decoded_file))
                smiles_list = [row[0].strip() for row in reader if row and row[0].strip()]
            except Exception as e:
                return Response({"error": f"Failed to parse CSV: {e}"}, status=400)
        elif smiles_input:
            if isinstance(smiles_input, str):
                smiles_list = [s.strip() for s in smiles_input.split(",") if s.strip()]
            elif isinstance(smiles_input, list):
                smiles_list = [s.strip() for s in smiles_input if isinstance(s, str)]
        else:
            return Response(
                {"error": "Provide either 'smiles' or 'file' (CSV)."},
                status=400,
            )

        # 🔹 Hapus duplikat
        smiles_list = list(dict.fromkeys(smiles_list))

        if not smiles_list:
            return Response({"error": "No valid SMILES provided."}, status=400)

        # 🔹 Pisahkan valid dan invalid SMILES
        valid_smiles, invalid_smiles = [], []
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            if mol:
                valid_smiles.append(s)
            else:
                invalid_smiles.append(s)

        predictions_dict = {s: None for s in smiles_list}
        ml_model = get_object_or_404(
            MLModel, descriptor=model_descriptor, method=model_method, is_active=True
        )

        try:
            # 🔹 Prediksi hanya untuk valid SMILES
            if valid_smiles:
                response = requests.post(
                    # "http://localhost:8080/api/v1/predict/",
                    settings.ML_MODEL_URL + "/api/v1/predict/",
                    json={
                        "smiles": valid_smiles,
                        "model_method": model_method,
                        "model_descriptor": model_descriptor,
                    },
                )

                if response.status_code != 200:
                    return Response(
                        {"error": "ML API error", "details": response.text}, status=400
                    )

                predictions_response = response.json()

                # Kalau ML API pakai format {"results": [...]}
                if "results" in predictions_response:
                    predictions = predictions_response["results"]
                else:
                    predictions = predictions_response  # fallback kalau format lama

                for item in predictions:
                    smiles = item.get("smiles")
                    pred_value = item.get("prediction")
                    predictions_dict[smiles] = pred_value

            # 🔹 Buat objek Prediction utama
            prediction = Prediction.objects.create(
                user=user,
                ml_model=ml_model,
                input_source_type="csv" if csv_file else "text",
                completed_at=timezone.now(),
            )

            # 🔹 Ambil compound dari DB (untuk yang valid)
            existing = Compound.objects.filter(smiles__in=valid_smiles)
            compound_map = {c.smiles: c for c in existing}
            compounds_to_fetch = [s for s in valid_smiles if s not in compound_map]

            # 🔹 Fetch compound baru dari PubChem
            fetched = {}
            if compounds_to_fetch:
                with ThreadPoolExecutor() as executor:
                    futures = {
                        executor.submit(self.fetch_pubchem_data, s): s
                        for s in compounds_to_fetch
                    }
                    for future in as_completed(futures):
                        s = futures[future]
                        try:
                            data = future.result()
                            fetched[s] = Compound.objects.create(smiles=s, **data)
                        except Exception as e:
                            print(f"⚠️ PubChem fetch failed for {s}: {e}")
                            fetched[s] = None

            all_compounds = {**compound_map, **fetched}

            # 🔹 Siapkan list untuk bulk_create
            to_create = []
            response_data = []

            for smiles in smiles_list:
                ic50 = predictions_dict.get(smiles)
                mol = Chem.MolFromSmiles(smiles)
                compound = all_compounds.get(smiles) if smiles in valid_smiles else None

                if mol and ic50 is not None:
                    # Hitung LELP dan kategori
                    heavy_atoms = mol.GetNumHeavyAtoms()
                    clogP = Crippen.MolLogP(mol)
                    le = (1.37 * ic50) / heavy_atoms if heavy_atoms else None
                    lelp = clogP / le if le else None
                    # Determine category based on IC50 (µM) thresholds:
                    # IC50 ≤ 20 -> high potential, 20 < IC50 ≤ 100 -> moderate, IC50 > 100 -> low potential
                    if ic50 is not None:
                        try:
                            ic50_val = float(ic50)
                        except Exception:
                            ic50_val = None

                        if ic50_val is not None:
                            if ic50_val >= 6.0:
                                category = "high potential"
                            elif ic50_val >= 5.0:
                                category = "moderate"
                            else:
                                category = "low potential"
                        else:
                            category = None
                    else:
                        category = None
                else:
                    ic50 = None
                    lelp = None
                    category = None

                # 🔹 Simpan semuanya, termasuk invalid
                pred_comp = PredictionCompound(
                    prediction=prediction,
                    compound=compound,
                    ic50=ic50,
                    lelp=lelp,
                    category=category,
                )
                to_create.append(pred_comp)

                response_data.append({
                    "smiles": smiles,
                    "pic50": ic50,
                    "lelp": lelp,
                    "category": category,
                    "compound": {
                        "id": compound.id if compound else None,
                        "smiles": smiles,
                        "iupac_name": getattr(compound, "iupac_name", None),
                        "cid": getattr(compound, "cid", None),
                        "description": getattr(compound, "description", None),
                        "molecular_formula": getattr(compound, "molecular_formula", None),
                        "molecular_weight": getattr(compound, "molecular_weight", None),
                        "synonyms": getattr(compound, "synonyms", None),
                        "inchi": getattr(compound, "inchi", None),
                        "inchikey": getattr(compound, "inchikey", None),
                        "structure_image": getattr(compound, "structure_image", None),
                    } if compound else None,
                })

            # 🔹 Simpan semua hasil ke DB
            PredictionCompound.objects.bulk_create(to_create)

            return Response({
                "status": "success",
                "message": "Prediction completed.",
                "valid_count": len(valid_smiles),
                "invalid_count": len(invalid_smiles),
                "invalid_smiles": invalid_smiles,
                "data": response_data,
            }, status=200)

        except Exception as e:
            return Response({"error": str(e)}, status=400)


    def fetch_pubchem_data(self, smiles):
        """Ambil data compound dari PubChem."""
        data = {
            "cid": None, "molecular_formula": None, "molecular_weight": None,
            "iupac_name": None, "inchi": None, "inchikey": None,
            "description": None, "synonyms": None, "structure_image": None,
        }
        try:
            compounds = pcp.get_compounds(smiles, "smiles")
            if compounds:
                c = compounds[0]
                data.update({
                    "cid": c.cid,
                    "molecular_formula": c.molecular_formula,
                    "molecular_weight": c.molecular_weight,
                    "iupac_name": c.iupac_name,
                    "inchi": c.inchi,
                    "inchikey": c.inchikey,
                    "synonyms": ", ".join(c.synonyms) if c.synonyms else None,
                    "structure_image": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={c.cid}&t=l" if c.cid else None,
                })
                desc = self.fetch_pubchem_description(c.cid)
                if desc:
                    data["description"] = desc
        except Exception as e:
            print(f"⚠️ PubChem fetch error for {smiles}: {e}")
        return data

    def fetch_pubchem_description(self, cid):
        try:
            res = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON",
                timeout=5,
            )
            if res.status_code == 200:
                info = res.json().get("InformationList", {}).get("Information", [])
                for i in info:
                    if "Description" in i:
                        return i["Description"]
        except Exception:
            pass
        return None

@extend_schema(**prediction_download_schema)
class PredictionDownloadView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request, id, *args, **kwargs):
        try:
            prediction = Prediction.objects.get(id=id)
        except Prediction.DoesNotExist:
            return Response(
                {
                    "status": "error",
                    "message": "Prediction not found.",
                    "data": None
                },
                status=status.HTTP_404_NOT_FOUND
            )

        if prediction.user != request.user and not getattr(request.user, "is_admin", False):
            return Response(
                {
                    "status": "error",
                    "message": "You do not have permission to download this prediction.",
                    "data": None
                },
                status=status.HTTP_403_FORBIDDEN
            )

        response = HttpResponse(content_type="text/csv")
        response["Content-Disposition"] = f'attachment; filename="prediction_{prediction.id}.csv"'

        writer = csv.writer(response)
        writer.writerow([
            "SMILES",
            "IUPAC_Name",
            "Molecular_Formula",
            "Molecular_Weight",
            "Predicted_IC50",
            "Predicted_Category",
            "PubChem_CID",
        ])

        prediction_results = prediction.prediction_compounds.select_related("compound").all()
        for result in prediction_results:
            writer.writerow([
                result.compound.smiles,
                result.compound.iupac_name,
                result.compound.molecular_formula,
                result.compound.molecular_weight,
                result.ic50,
                result.category,
                result.compound.cid,
            ])

        return response