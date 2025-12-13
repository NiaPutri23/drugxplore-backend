from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from rest_framework.decorators import action
from rest_framework.parsers import MultiPartParser, JSONParser
from drf_spectacular.utils import extend_schema_view, extend_schema
from django.db import transaction
import csv, math
from rdkit import Chem
from rdkit.Chem import Crippen

from api.models import PredictionCompound, LiteratureCompound, Compound
from .serializers import PredictionCompoundSerializer, LiteratureCompoundSerializer
from .permissions import IsAdminOrReadOnly
from api.v1.utils.pubchem_utils import get_pubchem_data

from .schemas import (
    prediction_compound_list_responses,
    prediction_compound_retrieve_responses,
    prediction_compound_destroy_responses,
    prediction_compound_lib_responses,
)


# ==============================================================
# 🔹 PREDICTION COMPOUND VIEWSET
# ==============================================================
@extend_schema_view(
    list=extend_schema(description="List all prediction compounds", responses=prediction_compound_list_responses),
    retrieve=extend_schema(description="Retrieve a specific prediction compound", responses=prediction_compound_retrieve_responses),
    destroy=extend_schema(description="Delete a specific prediction compound", responses=prediction_compound_destroy_responses),
)
class PredictionCompoundViewSet(viewsets.ModelViewSet):
    serializer_class = PredictionCompoundSerializer
    permission_classes = [IsAuthenticated]
    http_method_names = ['get', 'delete', 'head', 'options']

    def get_queryset(self):
        qs = PredictionCompound.objects.select_related('prediction', 'prediction__user', 'compound')
        if self.request.user.role != 'admin':
            qs = qs.filter(prediction__user=self.request.user)
        return qs.order_by('compound').distinct('compound')

    def list(self, request, *args, **kwargs):
        serializer = self.get_serializer(self.get_queryset(), many=True)
        return Response({
            "status": "success",
            "message": "Prediction compounds retrieved successfully.",
            "data": serializer.data,
        })

    def retrieve(self, request, *args, **kwargs):
        serializer = self.get_serializer(self.get_object())
        return Response({
            "status": "success",
            "message": "Prediction compound retrieved successfully.",
            "data": serializer.data,
        })

    def destroy(self, request, *args, **kwargs):
        self.get_object().delete()
        return Response({
            "status": "success",
            "message": "Prediction compound deleted successfully.",
        })


# ==============================================================
# 🔹 LITERATURE COMPOUND VIEWSET
# ==============================================================
@extend_schema(
    description="Retrieve or manage literature compound data.",
    responses=prediction_compound_lib_responses,
)
class LiteratureCompoundViewSet(viewsets.ModelViewSet):
    queryset = LiteratureCompound.objects.select_related('compound').all()
    serializer_class = LiteratureCompoundSerializer
    permission_classes = [IsAdminOrReadOnly]
    parser_classes = [MultiPartParser, JSONParser]

    # ==============================================================
    # 🔹 BULK CREATE MANUAL
    # ==============================================================
    @action(detail=False, methods=['post'], permission_classes=[IsAdminOrReadOnly])
    @transaction.atomic
    def bulk_create(self, request):
        """
        Admin dapat menambahkan beberapa compound sekaligus (manual input)
        Format JSON: [{ "smiles": "...", "ic50": ... }, ...]
        """
        compounds = request.data
        if not isinstance(compounds, list):
            return Response({"error": "Expected a list of compounds."}, status=400)

        created, skipped = 0, 0

        for item in compounds:
            smiles = item.get("smiles")
            ic50 = item.get("ic50")

            if not smiles or not ic50:
                skipped += 1
                continue

            try:
                ic50 = float(ic50)
                pic50 = -math.log10(ic50 * 1e-6)  # µM → M
            except ValueError:
                skipped += 1
                continue

            # Cek existing
            compound, _ = Compound.objects.get_or_create(smiles=smiles)
            if LiteratureCompound.objects.filter(compound=compound).exists():
                skipped += 1
                continue

            # Hitung LELP
            lelp, category = None, None
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    heavy_atoms = mol.GetNumHeavyAtoms()
                    clogP = Crippen.MolLogP(mol)
                    if heavy_atoms > 0:
                        le = (1.37 * ic50) / heavy_atoms
                        lelp = clogP / le if le else None
                        if lelp is not None:
                            category = (
                                "low potential" if lelp > 20 else
                                "moderate" if lelp >= 10 else
                                "potential"
                            )
            except Exception as e:
                print(f"⚠️ LELP calc failed for {smiles}: {e}")

            LiteratureCompound.objects.create(
                compound=compound,
                ic50=ic50,
                ic50val=pic50,
                lelp=lelp,
                category=category,
            )
            created += 1

        return Response({
            "status": "success",
            "message": f"{created} compounds added, {skipped} skipped (duplicate or invalid)."
        }, status=status.HTTP_201_CREATED)

    # ==============================================================
    # 📂 Upload CSV (Hanya smiles + ic50)
    # ==============================================================
    @action(detail=False, methods=['post'], permission_classes=[IsAdminOrReadOnly])
    @transaction.atomic
    def upload_csv(self, request):
        """
        Upload CSV dengan kolom:
        - smiles
        - ic50 (µM)
        """
        csv_file = request.FILES.get("file")
        if not csv_file:
            return Response({"error": "No CSV file uploaded."}, status=400)

        decoded = csv_file.read().decode("utf-8").splitlines()
        reader = csv.DictReader(decoded)

        created, skipped = 0, 0

        for row in reader:
            smiles = row.get("smiles")
            ic50_val = row.get("ic50")
            if not smiles or not ic50_val:
                skipped += 1
                continue

            try:
                ic50 = float(ic50_val)
                pic50 = -math.log10(ic50 * 1e-6)
            except ValueError:
                skipped += 1
                continue

            compound, _ = Compound.objects.get_or_create(smiles=smiles)
            if LiteratureCompound.objects.filter(compound=compound).exists():
                skipped += 1
                continue

            # Hitung LELP
            lelp, category = None, None
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    heavy_atoms = mol.GetNumHeavyAtoms()
                    clogP = Crippen.MolLogP(mol)
                    if heavy_atoms > 0:
                        le = (1.37 * ic50) / heavy_atoms
                        lelp = clogP / le if le else None
                        if lelp is not None:
                            category = (
                                "low potential" if lelp > 20 else
                                "moderate" if lelp >= 10 else
                                "high potential"
                            )
            except Exception as e:
                print(f"⚠️ LELP calc failed for {smiles}: {e}")

            LiteratureCompound.objects.create(
                compound=compound,
                ic50=ic50,
                ic50val=pic50,
                lelp=lelp,
                category=category,
            )
            created += 1

        return Response({
            "status": "success",
            "message": f"{created} records added, {skipped} skipped."
        }, status=status.HTTP_201_CREATED)