import importlib.util
import sys


def run_setup():
    modules = ['datetime', 'itertools', 'matplotlib', 'numpy', 'os', 'pandas', 'PIL', 'random',
               'rdkit', 'tkinter', 'tqdm', 'webbrowser']
    for mod in modules:
        if mod in sys.modules:
            print(f"{mod!r} already in sys.modules")
        elif (spec := importlib.util.find_spec(mod)) is not None:
            module = importlib.util.module_from_spec(spec)
            sys.modules[mod] = module
            spec.loader.exec_module(module)
            print(f"{mod!r} has been imported")
        else:
            print(f"can't find the {mod!r} module")