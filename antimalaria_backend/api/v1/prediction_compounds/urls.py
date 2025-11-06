from django.urls import path, include
from rest_framework.routers import DefaultRouter
from .views import PredictionCompoundViewSet, LiteratureCompoundViewSet

router = DefaultRouter()
# router.register(r'lib', LiteratureCompoundViewSet, basename='literature_compounds')
# router.register(r'', PredictionCompoundViewSet, basename='prediction_compounds')

router.register(r'compounds', PredictionCompoundViewSet, basename='prediction_compounds')
router.register(r'lib', LiteratureCompoundViewSet, basename='literature_compounds')


urlpatterns = [
    path('', include(router.urls)),
]