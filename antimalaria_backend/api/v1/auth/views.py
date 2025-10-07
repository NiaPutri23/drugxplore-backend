from rest_framework_simplejwt.views import TokenObtainPairView
from .serializers import RegisterSerializer, LoginSerializer, ForgotPasswordSerializer, ResetPasswordSerializer
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.permissions import AllowAny
from drf_spectacular.utils import extend_schema
from .schemas import register_responses, register_examples, login_responses, login_examples, logout_responses
from django.contrib.auth.tokens import default_token_generator
from django.utils.http import urlsafe_base64_encode, urlsafe_base64_decode
from django.utils.encoding import force_bytes, force_str
from django.core.mail import send_mail
from django.conf import settings
from django.contrib.auth import get_user_model
User = get_user_model()

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
@extend_schema(
    request=LoginSerializer,
    responses=login_responses,
    examples=login_examples
)
class LoginView(TokenObtainPairView):
    serializer_class = LoginSerializer

    def post(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)

        data = serializer.validated_data

        return Response({
            "status": "success",
            "message": "Login successful",
            "data": data
        }, status=status.HTTP_200_OK)

@extend_schema(
    description="Logout endpoint to invalidate the refresh token.",
    responses=logout_responses
)
class LogoutView(APIView):
    def post(self, request, *args, **kwargs):
        return Response({
            "status": "success",
            "message": "Logout successful"
        }, status=status.HTTP_200_OK)

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