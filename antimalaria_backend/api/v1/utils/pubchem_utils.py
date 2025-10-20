from rdkit import Chem
from rdkit.Chem import inchi
import requests

def get_pubchem_data(smiles):
    """Cari data compound dari PubChem pakai SMILES → fallback ke InChIKey, termasuk synonyms dan 2D image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Invalid SMILES: {smiles}")
        return None

    canonical = Chem.MolToSmiles(mol, canonical=True)
    try:
        inchi_str = inchi.MolToInchi(mol)
        inchikey_str = inchi.InchiToInchiKey(inchi_str)
    except:
        inchi_str, inchikey_str = None, None

    props = "IUPACName,MolecularFormula,MolecularWeight,InChI,InChIKey,CanonicalSMILES"
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    pubchem_data = None

    # 🔹 Cari pakai SMILES
    try:
        res = requests.get(f"{base}/compound/smiles/{canonical}/property/{props}/JSON", timeout=10)
        data = res.json()
        if "PropertyTable" in data:
            pubchem_data = data["PropertyTable"]["Properties"][0]
            print(f"✅ Found by SMILES: {canonical}")
    except Exception as e:
        print(f"⚠️ PubChem SMILES query failed: {e}")

    # 🔹 Kalau gagal, coba pakai InChIKey
    if not pubchem_data and inchikey_str:
        try:
            res = requests.get(f"{base}/compound/inchikey/{inchikey_str}/property/{props}/JSON", timeout=10)
            data = res.json()
            if "PropertyTable" in data:
                pubchem_data = data["PropertyTable"]["Properties"][0]
                print(f"✅ Found by InChIKey: {inchikey_str}")
        except Exception as e:
            print(f"⚠️ PubChem InChIKey query failed: {e}")

    if pubchem_data:
        # 🔹 Ambil Synonyms
        try:
            res = requests.get(f"{base}/compound/cid/{pubchem_data['CID']}/synonyms/JSON", timeout=10)
            syn_data = res.json()
            if "InformationList" in syn_data:
                synonyms = syn_data["InformationList"]["Information"][0].get("Synonym", [])
                pubchem_data["Synonyms"] = "; ".join(synonyms)
        except Exception as e:
            pubchem_data["Synonyms"] = None
            print(f"⚠️ Failed to fetch synonyms: {e}")

        # 🔹 Ambil 2D Structure Image URL
        pubchem_data["StructureImage"] = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={pubchem_data['CID']}&t=l"

    return pubchem_data
