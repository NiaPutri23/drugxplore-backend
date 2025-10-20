import numpy as np
import pickle
import json
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from concurrent.futures import ThreadPoolExecutor
from django.conf import settings
from pathlib import Path

MODEL_DIR = settings.ML_MODEL_DIR
ECFP_JSON_FILE_PATH = Path(__file__).resolve().parent / "ecfp_features.json"
PUBCHEMFP_JSON_FILE_PATH = Path(__file__).resolve().parent / "pubchemfp_features.json"
MODELS = {}
FEATURES = {}

def load_all_models():
    MODEL_DIR.mkdir(parents=True, exist_ok=True) 

    for model_path in MODEL_DIR.iterdir():
        if not model_path.is_file():
            continue 

        try:
            if model_path.suffix == ".pkl":
                with model_path.open("rb") as f:
                    model = pickle.load(f)
            else:
                continue

            model_name = model_path.name
            MODELS[model_name] = model

        except Exception as e:
            print(f"Failed to load model {model_path.name}: {e}")

def load_model(model_name):
    model_path = MODEL_DIR / model_name
    if not model_path.is_file():
        raise FileNotFoundError(f"Model file '{model_name}' not found in '{MODEL_DIR}'.")

    try:
        if model_path.suffix == ".pkl":
            with model_path.open("rb") as f:
                model = pickle.load(f)
        else:
            raise ValueError("Unsupported model format.")
        
        model_name = model_path.name
        MODELS[model_name] = model

        print(f"Model '{model_name}' loaded successfully.")

    except Exception as e:
        print(f"Failed to load model {model_path.name}: {e}")
        return None

def load_features():
    try:
        with open(ECFP_JSON_FILE_PATH, "r") as f:
            ecfp_features = json.load(f)
            FEATURES["ecfp"] = ecfp_features
            print("Jumlah fitur ECFP:", len(FEATURES.get("ecfp", [])))
    except Exception as e:
        print(f"Failed to load features from {ECFP_JSON_FILE_PATH}: {e}")

    try:
        with open(PUBCHEMFP_JSON_FILE_PATH, "r") as f:
            pubchemfp_features = json.load(f)
            FEATURES["pubchemfp"] = pubchemfp_features
            print("Jumlah fitur PUBCHEMFP:", len(FEATURES.get("pubchemfp", [])))
    except Exception as e:
        print(f"Failed to load features from {PUBCHEMFP_JSON_FILE_PATH}: {e}")

FEATURIZER_MAP = {
    "ECFP": lambda smiles, model=None: smiles_to_ecfp(smiles, tuple(FEATURES.get("ecfp", [])), model=model),
    "PUBCHEMFP": lambda smiles, model=None: smiles_to_pubchemfp(smiles, tuple(FEATURES.get("pubchemfp", [])), model=model)
}

def smiles_to_ecfp(smiles, selected_features, model=None, radius=3, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None

    fp_generator = GetMorganGenerator(radius, fpSize=n_bits)
    fp = fp_generator.GetFingerprint(mol)
    array = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, array)

    indices = np.array(selected_features).flatten().astype(int)
    new_array = array[indices]

    # pad atau truncate sesuai model
    if model:
        expected_features = model.n_features_in_
        if len(new_array) < expected_features:
            # pad 0
            new_array = np.pad(new_array, (0, expected_features - len(new_array)), 'constant')
        elif len(new_array) > expected_features:
            # truncate
            new_array = new_array[:expected_features]

    return new_array

def smiles_to_pubchemfp(smiles, selected_features, model=None, radius=2, n_bits=881):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None

    fp_generator = GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = fp_generator.GetFingerprint(mol)
    array = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, array)

    indices = np.array(selected_features).flatten().astype(int)
    new_array = array[indices]

    # pad atau truncate sesuai model
    if model:
        expected_features = model.n_features_in_
        if len(new_array) < expected_features:
            new_array = np.pad(new_array, (0, expected_features - len(new_array)), 'constant')
        elif len(new_array) > expected_features:
            new_array = new_array[:expected_features]

    return new_array

def predict_batch_ic50(smiles_list, model_name, model_descriptor):
    model = MODELS.get(model_name)
    if model is None:
        raise ValueError(f"Model '{model_name}' not found or failed to load.")

    featurizer = FEATURIZER_MAP.get(model_descriptor)
    if featurizer is None:
        raise ValueError("Unsupported model descriptor")

    results = []
    for smiles in smiles_list:
        fp = featurizer(smiles, model=model)
        if fp is None:
            results.append({
                "smiles": smiles,
                "prediction": None,
                "status": "invalid"
            })
            continue
        results.append({
            "smiles": smiles,
            "features": fp,
            "status": "valid"
        })

    valid_entries = [r for r in results if r["status"] == "valid"]
    if valid_entries:
        try:
            fp_array = np.vstack([r["features"] for r in valid_entries])
            preds = model.predict(fp_array)
            for r, p in zip(valid_entries, preds):
                r["prediction"] = float(p)
        except Exception as e:
            for r in valid_entries:
                r["prediction"] = None
                r["status"] = f"error: {e}"
    else:
        return {
            "error": "No valid SMILES provided",
            "results": results
        }

    # buang fitur dari hasil akhir biar response bersih
    for r in results:
        r.pop("features", None)

    return {"results": results}