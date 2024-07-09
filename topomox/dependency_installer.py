# topomox/dependency_installer.py
import importlib
import subprocess
import sys

def install_and_import(package, import_name=None):
    import_name = import_name or package
    try:
        importlib.import_module(import_name)
    except ImportError:
        print(f"Installing {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        importlib.import_module(import_name)