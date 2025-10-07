from rest_framework_simplejwt.serializers import TokenObtainPairSerializer
from rest_framework import serializers
from django.contrib.auth import get_user_model
User = get_user_model()
from django.contrib.auth.tokens import default_token_generator
from django.utils.http import urlsafe_base64_encode, urlsafe_base64_decode
from django.utils.encoding import force_bytes, force_str


class RegisterSerializer(serializers.ModelSerializer):
    email = serializers.EmailField(required=True)
    password = serializers.CharField(write_only=True, required=True)
    password2 = serializers.CharField(write_only=True, required=True)

    class Meta:
        model = User
        fields = ['username', 'email', 'password', 'password2']
        read_only_fields = ['role']

    def validate(self, data):
        if data['password'] != data['password2']:
            raise serializers.ValidationError("Passwords must match.")
    
        if User.objects.filter(email=data['email']).exists():
            raise serializers.ValidationError({"email": "A user with that email already exists."})
        if User.objects.filter(username=data['username']).exists():
            raise serializers.ValidationError({"username": "A user with that username already exists."})
        return data

    def create(self, validated_data):
        validated_data.pop('password2')
        role = validated_data.pop('role', 'user')
        user = User.objects.create_user(**validated_data)
        user.role = role
        user.save(update_fields=['role'])
        return user

class LoginSerializer(TokenObtainPairSerializer):
    def validate(self, attrs):
        data = super().validate(attrs)
        data.update({
            "username": self.user.username,
            "email": self.user.email,
            "id": str(self.user.id),
            "role": self.user.role
        })
        return data

class ForgotPasswordSerializer(serializers.Serializer):
    email = serializers.EmailField()

class ResetPasswordSerializer(serializers.Serializer):
    password = serializers.CharField(write_only=True, required=True)
    password2 = serializers.CharField(write_only=True, required=True)

    def validate(self, data):
        if data['password'] != data['password2']:
            raise serializers.ValidationError("Passwords do not match.")
        return data