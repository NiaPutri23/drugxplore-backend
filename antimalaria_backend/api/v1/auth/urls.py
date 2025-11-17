from django.urls import path
from .views import RegisterView, LoginView, MeView,LogoutView, ForgotPasswordView, ResetPasswordView, CSRFTokenView

urlpatterns = [
    path('get-csrf-token/', CSRFTokenView.as_view(), name='csrf-token'),
    path('register/', RegisterView.as_view()),
    path('login/', LoginView.as_view()),
    path('me/', MeView.as_view()),
    path('logout/', LogoutView.as_view()),
    path('forgot-password/', ForgotPasswordView.as_view(), name='forgot-password'),
    path('reset-password/<uidb64>/<token>/', ResetPasswordView.as_view(), name='reset-password'),
]