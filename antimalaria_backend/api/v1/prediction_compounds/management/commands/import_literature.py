import csv
from rdkit import Chem
from rdkit.Chem import Crippen
from django.db import transaction
from django.core.management.base import BaseCommand
from api.models import Compound, LiteratureCompound
from api.v1.utils.pubchem_utils import get_pubchem_data


class Command(BaseCommand):
    help = "Import literature compounds from CSV and enrich data from PubChem"

    def add_arguments(self, parser):
        parser.add_argument('csv_path', type=str, help='Path to your CSV file')

    @transaction.atomic
    def handle(self, *args, **options):
        csv_path = options['csv_path']

        # 🔹 Hanya hapus LiteratureCompound, bukan Compound
        self.stdout.write(self.style.WARNING("⚠️ Menghapus semua data LiteratureCompound..."))
        LiteratureCompound.objects.all().delete()
        Compound.objects.all().delete()
        self.stdout.write(self.style.SUCCESS("✅ Semua data LiteratureCompound telah dihapus."))

        with open(csv_path, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                smiles = row.get('canonical_smiles')
                ic50 = float(row.get('standard_value') or 0)
                pic50 = float(row.get('pchembl_value') or 0)

                if not smiles:
                    continue

                # 🔹 Ambil atau buat Compound
                compound, created = Compound.objects.get_or_create(smiles=smiles)

                # 🔹 Kalau belum punya data lengkap, fetch dari PubChem
                if not compound.iupac_name:
                    pubchem_data = get_pubchem_data(smiles)
                    if pubchem_data:
                        compound.iupac_name = pubchem_data.get('IUPACName')
                        compound.cid = pubchem_data.get('CID')
                        compound.molecular_formula = pubchem_data.get('MolecularFormula')
                        compound.molecular_weight = pubchem_data.get('MolecularWeight')
                        compound.inchi = pubchem_data.get('InChI')
                        compound.inchikey = pubchem_data.get('InChIKey')
                        compound.synonyms = pubchem_data.get('Synonyms')
                        compound.structure_image = pubchem_data.get('StructureImage')
                        compound.save()

                # 🔹 Hitung LELP
                lelp, category = None, None
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol and ic50:
                        heavy_atoms = mol.GetNumHeavyAtoms()
                        clogP = Crippen.MolLogP(mol)
                        if heavy_atoms > 0:
                            le = (1.37 * ic50) / heavy_atoms
                            lelp = clogP / le if le else None
                            # Determine category based on IC50 (µM) thresholds:
                            # IC50 ≤ 20 -> high potential, 20 < IC50 ≤ 100 -> moderate, IC50 > 100 -> low potential
                            if ic50 and ic50 > 0:
                                if ic50 <= 20:
                                    category = "high potential"
                                elif ic50 <= 100:
                                    category = "moderate"
                                else:
                                    category = "low potential"
                            else:
                                category = None
                except Exception as e:
                    self.stdout.write(self.style.WARNING(f"⚠️ LELP calc failed for {smiles}: {e}"))

                # 🔹 Simpan LiteratureCompound
                LiteratureCompound.objects.create(
                    compound=compound,
                    ic50=ic50,
                    ic50val=pic50,
                    lelp=lelp,
                    category=category,
                )

                status = "🆕" if created else "🔁"
                # Include ic50 (µM) and ic50val (pIC50) in the import log
                self.stdout.write(self.style.SUCCESS(
                    f"{status} {compound.iupac_name or smiles} | IC50={ic50} µM | IC50val={pic50} | LELP={lelp} | Category={category}"
                ))

        self.stdout.write(self.style.SUCCESS("🎉 Import selesai tanpa hapus Compound!"))
