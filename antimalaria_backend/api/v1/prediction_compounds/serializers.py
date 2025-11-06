from rest_framework import serializers
from api.models import Compound, Prediction, PredictionCompound, MLModel, LiteratureCompound

class MLModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = MLModel
        fields = '__all__'

class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        exclude = ['created_at']

class PredictionSerializer(serializers.ModelSerializer):
    ml_model = MLModelSerializer(read_only=True)
    class Meta:
        model = Prediction
        fields = [
            "id",
            "user",
            "ml_model",
            "input_source_type",
            "created_at"
        ]

class PredictionCompoundSerializer(serializers.ModelSerializer):
    compound = CompoundSerializer(read_only=True)
    prediction = PredictionSerializer(read_only=True)
    class Meta:
        model = PredictionCompound
        fields = [
            'id', 'ic50', 'lelp', 'category', 'compound', 'prediction'
        ]

class LiteratureCompoundSerializer(serializers.ModelSerializer):
    compound = CompoundSerializer(read_only=True)
    compound_id = serializers.PrimaryKeyRelatedField(
        queryset=Compound.objects.all(),
        source='compound',
        write_only=True,
        required=False
    )

    class Meta:
        model = LiteratureCompound
        fields = [
            'id', 'compound', 'compound_id', 'ic50', 'ic50val', 'lelp', 'category',
        ]