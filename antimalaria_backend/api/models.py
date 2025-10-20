import uuid
from django.db import models
from django.conf import settings
from django.contrib.auth.models import AbstractUser

class CustomUser(AbstractUser):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    
    role = models.CharField(max_length=50, null=True, blank=True, default='user')  
    
    def __str__(self):
        return self.username

class MLModel(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=255, null=True, blank=True) 
    method = models.CharField(max_length=255, null=True, blank=True) 
    descriptor= models.CharField(max_length=255, null=True, blank=True) 
    version = models.CharField(max_length=50)  
    file = models.FileField(upload_to='ml_models/', null=True, blank=True)
    is_active = models.BooleanField(default=False) 
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.name} v{self.version}"

class ModelUpdate(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    ml_model = models.ForeignKey(MLModel, related_name='model_updates', on_delete=models.CASCADE)
    user = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='model_updates', on_delete=models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"Update for {self.ml_model.name} by {self.user.username} on {self.created_at}"

class Prediction(models.Model):
    class Meta:
        ordering = ['-created_at']

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE, null=True, blank=True)
    ml_model = models.ForeignKey(MLModel, on_delete=models.PROTECT, null=True, blank=True)
    input_source_type = models.CharField(max_length=50, null=True, blank=True) 
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return f"Prediction Job {self.id} ({self.status})"

class Compound(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    iupac_name = models.CharField(max_length=255, null=True, blank=True)  
    cid = models.CharField(max_length=50, null=True, blank=True) 
    smiles = models.TextField(null=True, blank=True) 
    description = models.TextField(null=True, blank=True) 
    molecular_formula = models.CharField(max_length=255, null=True, blank=True)  
    molecular_weight = models.FloatField(null=True, blank=True)   
    synonyms = models.TextField(null=True, blank=True)
    inchi = models.TextField(null=True, blank=True)  
    inchikey = models.CharField(max_length=255, null=True, blank=True)  
    structure_image = models.URLField(null=True, blank=True)  
    created_at = models.DateTimeField(auto_now_add=True) 

    def __str__(self):
        return self.iupac_name or "Unnamed Compound"

class PredictionCompound(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    prediction = models.ForeignKey(Prediction, related_name='prediction_compounds', on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, related_name='prediction_compounds', on_delete=models.CASCADE)

    ic50 = models.FloatField(null=True, blank=True)
    lelp = models.FloatField(null=True, blank=True)
    category = models.CharField(null=True, blank=True)
    class Meta:
        unique_together = ('prediction', 'compound')

    def __str__(self):
        return f"Result for {self.compound.iupac_name} in Job {self.prediction.id}"
    
class LiteratureCompound(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    compound = models.ForeignKey(Compound, related_name='literature_data', on_delete=models.CASCADE)
    
    ic50 = models.FloatField(null=True, blank=True)
    ic50val = models.FloatField(null=True, blank=True)
    lelp = models.FloatField(null=True, blank=True)
    category = models.CharField(max_length=100, null=True, blank=True)

    def __str__(self):
        return f"Literature data for {self.compound.iupac_name}"