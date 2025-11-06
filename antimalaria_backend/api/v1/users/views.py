from rest_framework import status, viewsets
from .serializers import UserSerializer
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth import get_user_model
from rest_framework.exceptions import PermissionDenied
from django.contrib.auth.password_validation import validate_password
from django.core.exceptions import ValidationError
from drf_spectacular.utils import extend_schema, extend_schema_view, OpenApiResponse

User = get_user_model()


@extend_schema_view(
    list=extend_schema(
        description="Get a list of users (admin only). Regular users only see their own data.",
        responses={
            200: OpenApiResponse(
                description="List of users.",
                response={
                    "type": "object",
                    "properties": {
                        "status": {"type": "string"},
                        "message": {"type": "string"},
                        "data": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "uuid"},
                                "username": {"type": "string"},
                                "email": {"type": "string"},
                                "role": {"type": "string"},
                            }
                        },
                    }
                }
            ),
            403: OpenApiResponse(description="Forbidden: Not allowed."),
        }
    ),
    retrieve=extend_schema(
        description="Retrieve a user's details. Regular users can only access their own data.",
        responses={
            200: OpenApiResponse(
                description="User details.",
                response={
                    "type": "object",
                    "properties": {
                        "status": {"type": "string"},
                        "message": {"type": "string"},
                        "data": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "string"},
                                "username": {"type": "string"},
                                "email": {"type": "string"},
                                "role": {"type": "string"},
                            }
                        },
                    }
                }
            ),
            403: OpenApiResponse(description="Forbidden: Not allowed."),
            404: OpenApiResponse(description="User not found."),
        }
    ),
    destroy=extend_schema(
        description="Delete a user (admin only).",
        responses={
            204: OpenApiResponse(
                description="User deleted successfully.",
                response={
                    "type": "object",
                    "properties": {
                        "status": {"type": "string"},
                        "message": {"type": "string"},
                    }
                }
            ),
            403: OpenApiResponse(description="Forbidden: Not allowed."),
            404: OpenApiResponse(description="User not found.")
        }
    ),
    partial_update=extend_schema(
        description="Update the password for the current user. Requires both `password` and `password2` fields.",
        request={
            "application/json": {
                "type": "object",
                "properties": {
                    "password": {"type": "string"},
                    "password2": {"type": "string"},
                },
                "required": ["password", "password2"]
            }
        },
        responses={
            200: OpenApiResponse(
                description="Password has been changed.",
                response={
                    "type": "object",
                    "properties": {
                        "status": {"type": "string"},
                        "message": {"type": "string"},
                    }
                }
            ),
            400: OpenApiResponse(description="Validation error."),
            403: OpenApiResponse(description="Forbidden: Not allowed."),
        }
    )
)
class UserViewSet(viewsets.ModelViewSet):
    serializer_class = UserSerializer
    permission_classes = [IsAuthenticated]
    http_method_names = ['get', 'patch', 'delete', 'head', 'options']

    def get_queryset(self):
        if self.request.user.role != 'admin':
            return User.objects.filter(id=self.request.user.id)
        return User.objects.all()

    def get_object(self):
        obj = super().get_object()
        if self.request.user.role != 'admin' and obj.id != self.request.user.id:
            raise PermissionDenied("You do not have permission to access this user.")
        return obj
    
    def list(self, request, *args, **kwargs):
        queryset = self.filter_queryset(self.get_queryset())
        serializer = self.get_serializer(queryset, many=True)
        return Response({
            "status": "success",
            "message": "Users retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)
    
    def retrieve(self, request, *args, **kwargs):
        instance = self.get_object()
        serializer = self.get_serializer(instance)
        return Response({
            "status": "success",
            "message": "User retrieved successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)

    def destroy(self, request, *args, **kwargs):
        if request.user.role != 'admin':
            raise PermissionDenied("Only admin can delete users.")
        instance = self.get_object()
        self.perform_destroy(instance)
        return Response({
            "status": "success",
            "message": "User deleted successfully."
        }, status=status.HTTP_204_NO_CONTENT)
    
def partial_update(self, request, *args, **kwargs):
    user = self.get_object()

    # 🚧 Batasi akses
    if request.user.role != 'admin' and user.id != request.user.id:
        raise PermissionDenied("You do not have permission to update this user.")

    password = request.data.get("password")
    password2 = request.data.get("password2")
    username = request.data.get("username")

    # 🔹 Jika ubah password
    if password or password2:
        if not password or not password2:
            return Response(
                {"detail": "Both password and password2 are required."},
                status=status.HTTP_400_BAD_REQUEST
            )

        if password != password2:
            return Response(
                {"detail": "Passwords do not match."},
                status=status.HTTP_400_BAD_REQUEST
            )

        try:
            validate_password(password, user)
        except ValidationError as e:
            return Response({"detail": e.messages}, status=status.HTTP_400_BAD_REQUEST)

        user.set_password(password)
        user.save()
        return Response({
            "status": "success",
            "message": "Password has been changed."
        }, status=status.HTTP_200_OK)

    # 🔹 Jika ubah username
    if username:
        # Pastikan username unik
        if User.objects.filter(username=username).exclude(id=user.id).exists():
            return Response(
                {"detail": "This username is already taken."},
                status=status.HTTP_400_BAD_REQUEST
            )

        user.username = username
        user.save()

        serializer = self.get_serializer(user)
        return Response({
            "status": "success",
            "message": "Username updated successfully.",
            "data": serializer.data
        }, status=status.HTTP_200_OK)

    # 🔹 Tidak ada field valid yang dikirim
    return Response(
        {"detail": "No valid fields provided. You can only update username or password."},
        status=status.HTTP_400_BAD_REQUEST
    )
