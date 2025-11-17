from django.contrib.auth import get_user_model, authenticate, login, logout
from django.contrib.auth.tokens import default_token_generator
from django.utils.http import urlsafe_base64_encode, urlsafe_base64_decode
from django.utils.encoding import force_bytes, force_str
from django.core.mail import send_mail
from django.conf import settings
from django.middleware.csrf import get_token

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status, serializers
from rest_framework.permissions import AllowAny, IsAuthenticated

from .serializers import RegisterSerializer, ForgotPasswordSerializer, ResetPasswordSerializer
from .schemas import register_responses, register_examples, login_responses, login_examples, logout_responses
from drf_spectacular.utils import extend_schema

User = get_user_model()

# Serializer sederhana untuk login sesi
class SessionLoginSerializer(serializers.Serializer):
    username = serializers.CharField(required=True)
    password = serializers.CharField(required=True, write_only=True)

# View untuk frontend mengambil CSRF token saat pertama kali load
class CSRFTokenView(APIView):
    permission_classes = [AllowAny]

    def get(self, request, *args, **kwargs):
        # Memastikan cookie csrftoken diatur
        get_token(request) 
        return Response({"detail": "CSRF cookie set."})

# Register View
@extend_schema(
    description="User registration endpoint.",
    request=RegisterSerializer,
    responses=register_responses,
    examples=register_examples
)
class RegisterView(APIView):
    permission_classes = [AllowAny] 

    def post(self, request, *args, **kwargs):
        serializer = RegisterSerializer(data=request.data)
        if serializer.is_valid():
            user = serializer.save()
            return Response({
                "status": "success",
                "message": "User registered successfully.",
                "data": {
                    "id": str(user.id),
                    "username": user.username,
                    "email": user.email,
                    "role": user.role
                    }
                },
                status=status.HTTP_201_CREATED
            )
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

# Login view menggunakan session authentication
@extend_schema(
    request=SessionLoginSerializer, # Menggunakan serializer baru
    responses=login_responses,
    examples=login_examples
)
class LoginView(APIView):
    permission_classes = [AllowAny]

    def post(self, request, *args, **kwargs):
        serializer = SessionLoginSerializer(data=request.data)
        if not serializer.is_valid():
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        username = serializer.validated_data['username']
        password = serializer.validated_data['password']
        
        # Menggunakan authenticate bawaan Django
        user = authenticate(request, username=username, password=password)

        if user is not None:
            # Menggunakan login bawaan Django, ini akan SETEL cookie sessionid
            login(request, user)
            
            # Kita kirim data user kembali, mirip seperti RegisterView
            return Response({
                "status": "success",
                "message": "Login successful",
                "data": {
                    "id": str(user.id),
                    "username": user.username,
                    "email": user.email,
                    "role": user.role
                    }
            }, status=status.HTTP_200_OK)
        else:
            return Response({
                "status": "error",
                "message": "Invalid username or password"
            }, status=status.HTTP_401_UNAUTHORIZED)
        

class MeView(APIView):
    permission_classes = [IsAuthenticated]

    def get(self, request):
        user = request.user
        return Response({
            "id": user.id,
            "username": user.username,
            "email": user.email,
            "role": user.role
        })

# Logout view
@extend_schema(
    description="Logout endpoint untuk menghapus session cookie.",
    responses=logout_responses
)
class LogoutView(APIView):
    # Hanya user yang terautentikasi (punya cookie) yang bisa logout
    permission_classes = [IsAuthenticated] 

    def post(self, request, *args, **kwargs):
        # Menggunakan logout bawaan Django, ini akan HAPUS cookie sessionid
        logout(request)
        
        return Response({
            "status": "success",
            "message": "Logout successful"
        }, status=status.HTTP_200_OK)

# ForgotPasswordView
class ForgotPasswordView(APIView):
    permission_classes = [AllowAny]

    def post(self, request):
        serializer = ForgotPasswordSerializer(data=request.data)
        if serializer.is_valid():
            email = serializer.validated_data['email']
            try:
                user = User.objects.get(email=email)
            except User.DoesNotExist:
                return Response({"error": "User not found."}, status=status.HTTP_404_NOT_FOUND)
            uid = urlsafe_base64_encode(force_bytes(user.pk))
            token = default_token_generator.make_token(user)
            reset_url = f"{settings.FRONTEND_URL}/reset-password/{uid}/{token}/"
            send_mail(
                "Password Reset",
                f"Click the link to reset your password: {reset_url}",
                settings.DEFAULT_FROM_EMAIL,
                [email],
            )
            return Response({"message": "Password reset link sent. Check your email"}, status=status.HTTP_200_OK)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

# ResetPasswordView
class ResetPasswordView(APIView):
    permission_classes = [AllowAny]

    def post(self, request, uidb64, token):
        serializer = ResetPasswordSerializer(data=request.data)
        if serializer.is_valid():
            try:
                uid = force_str(urlsafe_base64_decode(uidb64))
                user = User.objects.get(pk=uid)
            except (TypeError, ValueError, OverflowError, User.DoesNotExist):
                return Response({"error": "Invalid link."}, status=status.HTTP_400_BAD_REQUEST)
            if not default_token_generator.check_token(user, token):
                return Response({"error": "Invalid or expired token."}, status=status.HTTP_400_BAD_REQUEST)
            user.set_password(serializer.validated_data['password'])
            user.save()
            return Response({"message": "Password has been reset."}, status=status.HTTP_200_OK)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)