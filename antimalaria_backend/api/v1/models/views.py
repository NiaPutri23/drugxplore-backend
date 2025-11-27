from rest_framework import status, viewsets
from antimalaria_backend import settings
from .serializers import MLModelSerializer
from .permissions import IsAdminUser
from api.models import MLModel, ModelUpdate
from django.shortcuts import get_object_or_404
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser, FormParser
from rest_framework.views import APIView
from django.conf import settings
# import environ
# env = environ.Env()
# environ.Env.read_env()
import requests
from drf_spectacular.utils import extend_schema, extend_schema_view
from .schemas import model_create_responses, model_list_responses, model_retrieve_responses, model_delete_responses, model_activate_responses

@extend_schema_view(
    list=extend_schema(
        description="Get a list of ML models (admin only).",
        responses=model_list_responses
    ),
    retrieve=extend_schema(
        description="Retrieve a specific ML model by ID (admin only).",
        responses=model_retrieve_responses
    ),
    create=extend_schema(
        description="Upload a new ML model (admin only).",
        responses=model_create_responses,
    ),
    destroy=extend_schema(
        description="Delete an ML model by ID (admin only).",
        responses=model_delete_responses
    )
)
class MLViewSet(viewsets.ModelViewSet):
    serializer_class = MLModelSerializer
    permission_classes = [IsAdminUser]
    parser_classes = [MultiPartParser, FormParser]
    http_method_names = ['get', 'post', 'delete', 'head', 'options']

    def get_queryset(self):
        return MLModel.objects.all().order_by('-created_at')
    
    def list(self, request, *args, **kwargs):
        queryset = self.get_queryset()
        serializer = self.get_serializer(queryset, many=True)
        return Response({
            "status": "success",
            "message": "Models retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)
    
    def retrieve(self, request, *args, **kwargs):
        instance = self.get_object()
        serializer = self.get_serializer(instance)
        return Response({
            "status": "success",
            "message": "Model retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)

    def create(self, request, *args, **kwargs):
        descriptor = request.data.get("model_descriptor")
        method = request.data.get("model_method")

        if not descriptor or not method:
            return Response({"error": "model_descriptor and model_method are required."}, status=status.HTTP_400_BAD_REQUEST)

        if not model_file:
            return Response({"error": "No file provided."}, status=status.HTTP_400_BAD_REQUEST)
        
        if not model_file.name.endswith('.pkl'):
            return Response({"error": "Invalid file type. Only .pkl files are allowed."}, status=status.HTTP_400_BAD_REQUEST)

        latest_model = MLModel.objects.filter(
            descriptor=descriptor,
            method=method
        ).order_by('-version').first()
        new_version = int(latest_model.version) + 1 if latest_model else 1

        model_file = request.FILES.get("file")
        backend_url = settings.ML_MODEL_URL


        model_file.seek(0)

        data = {
            "name": model_file.name,
            "method": method,
            "descriptor": descriptor,
            "version": str(new_version),
            "is_active": "false",
        }

        response = requests.post(backend_url, data=data, files={"file": model_file})

        if response.status_code != 201:
            return Response(
                {"error": "Failed to store model in backend_model service",
                "backend_response": response.text},
                status=status.HTTP_502_BAD_GATEWAY
            )

        return Response({
            "status": "success",
            "message": "Model uploaded successfully.",
            "data": response.json()
        }, status=status.HTTP_201_CREATED)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        backend_url = f"{settings.ML_MODEL_URL}{instance.id}/"

        response = requests.delete(backend_url)

        if response.status_code not in (200, 204):
            return Response(
                {"error": "Failed to delete model in backend_model service",
                "backend_response": response.text},
                status=status.HTTP_502_BAD_GATEWAY
            )
        instance.delete()
        return Response({
            "status": "success",
            "message": f"Model {instance.name} deleted successfully."
        }, status=status.HTTP_200_OK)


@extend_schema(
    description="Activate a specific ML model by ID (admin only). Deactivates other models with the same descriptor and method.",
    responses=model_activate_responses
)
class ActivateModelAPIView(APIView):
    permission_classes = [IsAdminUser]

    def post(self, request, *args, **kwargs):
        id = kwargs.get("id")
        model = get_object_or_404(MLModel, id=id)
        model.is_active = True
        model.save()
        MLModel.objects.filter(
            descriptor=model.descriptor,
            method=model.method
        ).exclude(id=model.id).update(is_active=False)
        return Response({
            "status": "success",
            "message": f"Model {model.name} activated successfully."
        }, status=status.HTTP_200_OK)