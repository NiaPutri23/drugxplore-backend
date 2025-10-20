from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated, AllowAny
from drf_spectacular.utils import extend_schema_view, extend_schema
from rest_framework.views import APIView
from .serializers import PredictionCompoundSerializer, LiteratureCompoundSerializer
from api.models import PredictionCompound, LiteratureCompound
from .permissions import IsAdminOrReadOnly
from .schemas import (
    prediction_compound_list_responses,
    prediction_compound_retrieve_responses,
    prediction_compound_destroy_responses,
    prediction_compound_lib_responses
)
@extend_schema_view(
    list=extend_schema(description="List all prediction compounds", responses=prediction_compound_list_responses),
    retrieve=extend_schema(description="Retrieve a specific prediction compound", responses=prediction_compound_retrieve_responses),
    destroy=extend_schema(description="Delete a specific prediction compound", responses=prediction_compound_destroy_responses)
)
class PredictionCompoundViewSet(viewsets.ModelViewSet):
    serializer_class = PredictionCompoundSerializer
    permission_classes = [IsAuthenticated]
    http_method_names = ['get', 'delete', 'head', 'options']

    def get_queryset(self):
        qs = PredictionCompound.objects.select_related('prediction', 'prediction__user', 'compound')
        if self.request.user.role != 'admin':
            qs = qs.filter(prediction__user=self.request.user)
        
        qs = qs.order_by('compound').distinct('compound')
        return qs

    def list(self, request, *args, **kwargs):
        queryset = self.get_queryset()
        serializer = self.get_serializer(queryset, many=True)
        return Response({
            "status": "success",
            "message": "Prediction compounds retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)

    def retrieve(self, request, *args, **kwargs):
        instance = self.get_object()
        serializer = self.get_serializer(instance)
        return Response({
            "status": "success",
            "message": "Prediction compound retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        instance.delete()
        return Response({
            "status": "success",
            "message": "Prediction compound deleted successfully."
        }, status=status.HTTP_200_OK)
    
@extend_schema(
    description="Retrieve a predefined library of prediction compounds.",
    responses=prediction_compound_lib_responses
)

# class PredictionCompoundLibView(APIView):
#     permission_classes = [AllowAny]

#     # dummy data
#     def get(self, request, *args, **kwargs):
#         data = [
#             # 1. Aspirin (Original)
#             {
#                 "ic50": 6.5,
#                 "lelp": 3.2,
#                 "category": "very strong",
#                 "compound": {
#                     "iupac_name": "acetylsalicylic acid",
#                     "cid": 2244,
#                     "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
#                     "description": "Aspirin is a medication used to reduce pain, fever, or inflammation.",
#                     "molecular_formula": "C9H8O4",
#                     "molecular_weight": 180.16,
#                     "synonyms": "2-(acetyloxy)benzoic acid, Aspirin",
#                     "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
#                     "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2244/PNG"
#                 }
#             },
#             # 2. Paracetamol (Original)
#             {
#                 "ic50": 5.1,
#                 "lelp": 5.6,
#                 "category": "strong",
#                 "compound": {
#                     "iupac_name": "paracetamol",
#                     "cid": 1983,
#                     "smiles": "CC(=O)NC1=CC=C(C=C1)O",
#                     "description": "Paracetamol, also known as acetaminophen, is a medication used to treat pain and fever.",
#                     "molecular_formula": "C8H9NO2",
#                     "molecular_weight": 151.16,
#                     "synonyms": "Acetaminophen, N-(4-hydroxyphenyl)acetamide",
#                     "inchi": "InChI=1S/C8H9NO2/c1-6(10)9-7-3-5-8(11)4-2-7/h2-5,11H,1H3,(H,9,10)",
#                     "inchikey": "RZVAJINKPMORJF-UHFFFAOYSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/1983/PNG"
#                 }
#             },
#             # 3. Ibuprofen
#             {
#                 "ic50": 4.5,
#                 "lelp": 2.1,
#                 "category": "moderate",
#                 "compound": {
#                     "iupac_name": "(RS)-2-(4-(2-methylpropyl)phenyl)propanoic acid",
#                     "cid": 3672,
#                     "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
#                     "description": "Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID) used for treating pain, fever, and inflammation.",
#                     "molecular_formula": "C13H18O2",
#                     "molecular_weight": 206.28,
#                     "synonyms": "Advil, Motrin",
#                     "inchi": "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
#                     "inchikey": "HEFNNWSXXWATRW-UHFFFAOYSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/3672/PNG"
#                 }
#             },
#             # 4. Metformin
#             {
#                 "ic50": 3.2,
#                 "lelp": 0.8,
#                 "category": "weak",
#                 "compound": {
#                     "iupac_name": "N,N-Dimethylimidodicarbonimidic diamide",
#                     "cid": 4091,
#                     "smiles": "CN(C)C(=N)N=C(N)N",
#                     "description": "Metformin is a first-line medication for the treatment of type 2 diabetes.",
#                     "molecular_formula": "C4H11N5",
#                     "molecular_weight": 129.16,
#                     "synonyms": "Glucophage",
#                     "inchi": "InChI=1S/C4H11N5/c1-9(2)4(7)8-3(5)6/h1-2H3,(H5,5,6,7,8)",
#                     "inchikey": "OETHQSOLBKYFAN-UHFFFAOYSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/4091/PNG"
#                 }
#             },
#             # 5. Caffeine
#             {
#                 "ic50": 3.7,
#                 "lelp": 0.9,
#                 "category": "weak",
#                 "compound": {
#                     "iupac_name": "1,3,7-Trimethylpurine-2,6-dione",
#                     "cid": 2519,
#                     "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
#                     "description": "Caffeine is a central nervous system stimulant of the methylxanthine class.",
#                     "molecular_formula": "C8H10N4O2",
#                     "molecular_weight": 194.19,
#                     "synonyms": "Guaranine, Theine",
#                     "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
#                     "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2519/PNG"
#                 }
#             },
#             # 6. Amoxicillin
#             {
#                 "ic50": 6.2,
#                 "lelp": 4.5,
#                 "category": "strong",
#                 "compound": {
#                     "iupac_name": "(2S,5R,6R)-6-[[(2R)-2-amino-2-(4-hydroxyphenyl)acetyl]amino]-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid",
#                     "cid": 33613,
#                     "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
#                     "description": "Amoxicillin is an antibiotic used to treat a number of bacterial infections.",
#                     "molecular_formula": "C16H19N3O5S",
#                     "molecular_weight": 365.4,
#                     "synonyms": "Amoxil, Moxatag",
#                     "inchi": "InChI=1S/C16H19N3O5S/c1-16(2)11(15(23)24)19-13(22)10(14(19)25-16)18-12(21)9(17)7-3-5-8(20)6-4-7/h3-6,9-11,14,20H,17H2,1-2H3,(H,18,21)(H,23,24)/t9-,10-,11+,14-/m1/s1",
#                     "inchikey": "LJQSCHXDOKSAOY-MNNPPOFNSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/33613/PNG"
#                 }
#             },
#             # 7. Lisinopril
#             {
#                 "ic50": 6.1,
#                 "lelp": 6.8,
#                 "category": "very strong",
#                 "compound": {
#                     "iupac_name": "(2S)-1-[(2S)-6-amino-2-[[(1S)-1-carboxy-3-phenylpropyl]amino]hexanoyl]pyrrolidine-2-carboxylic acid",
#                     "cid": 5362119,
#                     "smiles": "C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
#                     "description": "Lisinopril is a medication of the angiotensin-converting enzyme (ACE) inhibitor class used to treat high blood pressure, heart failure, and after heart attacks.",
#                     "molecular_formula": "C21H31N3O5",
#                     "molecular_weight": 405.49,
#                     "synonyms": "Zestril, Prinivil",
#                     "inchi": "InChI=1S/C21H31N3O5/c22-14-6-2-1-5-17(24-18(20(27)28)12-11-16-9-7-8-10-16)19(26)23-13-3-4-15(23)21(29)30/h7-10,15,17-18H,1-6,11-14,22H2,(H,24,26)(H,27,28)(H,29,30)/t15-,17-,18-/m0/s1",
#                     "inchikey": "ORRDEVOFBUBGPX-SSSRNJBGSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/5362119/PNG"
#                 }
#             },
#             # 8. Atorvastatin
#             {
#                 "ic50": 2.5,
#                 "lelp": 4.9,
#                 "category": "inactive",
#                 "compound": {
#                     "iupac_name": "(3R,5R)-7-[2-(4-fluorophenyl)-5-(propan-2-yl)-3-phenyl-4-(phenylcarbamoyl)-1H-pyrrol-1-yl]-3,5-dihydroxyheptanoic acid",
#                     "cid": 60823,
#                     "smiles": "CC(C)C1=CC(=C(N1CC[C@@H](C[C@@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C(=O)NC3=CC=CC=C3",
#                     "description": "Atorvastatin, sold under the brand name Lipitor among others, is a statin medication used to prevent cardiovascular disease and treat abnormal lipid levels.",
#                     "molecular_formula": "C33H35FN2O5",
#                     "molecular_weight": 558.64,
#                     "synonyms": "Lipitor",
#                     "inchi": "InChI=1S/C33H35FN2O5/c1-21(2)31-30(33(41)35-26-15-11-8-12-16-26)29(22-13-9-7-10-14-22)28(23-17-19-24(34)20-18-23)36(31)20-18-23)36(31)20-18-23)36(31)20-18-23)36(31)20-18-23",
#                     "inchikey": "XUKUURHRXDUEBC-KAYWLYCHSA-N",
#                     "structure_image": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/60823/PNG"
#                 }
#             }
#         ]
#         return Response({
#             "status": "success",
#             "message": "Prediction compound library retrieved successfully.",
#             "data": data
#         }, status=status.HTTP_200_OK)

class LiteratureCompoundViewSet(viewsets.ModelViewSet):
    queryset = LiteratureCompound.objects.all()
    serializer_class = LiteratureCompoundSerializer

    # 🔹 ini WAJIB ditulis bareng-bareng supaya GET bisa public
    authentication_classes = []  
    permission_classes = [IsAdminOrReadOnly]