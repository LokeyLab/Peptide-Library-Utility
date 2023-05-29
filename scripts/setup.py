import importlib.util
import sys


def run_module_check():
    modules = ['datetime', 'importlib', 'itertools', 'matplotlib', 'numpy', 'os',
               'pandas', 'PIL', 'pathlib2', 'pyshortcuts', 'random', 'rdkit', 'sys', 'tkinter',
               'tqdm', 'webbrowser']

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


def run_shortcut_creator():
    import os
    import platform
    import pathlib2
    from pyshortcuts import make_shortcut

    platform = platform.system().lower()
    print(f"\nPlatform is {platform}")

    if platform == "linux" or platform == "linux2":
        icon = pathlib2.Path(os.path.abspath("logo.ico"))

    elif platform == "darwin":
        icon = pathlib2.Path(os.path.abspath("logo.icns"))

    elif platform == "win32" or "win64":
        icon = pathlib2.Path(os.path.abspath("logo.ico"))

    main_script = pathlib2.Path(os.path.abspath("peptide_library_utility.py"))
    print(f"\n{main_script}")

    make_shortcut(str(main_script), name='Peptide Library Utility', icon=str(icon), desktop=True)

run_module_check()
# run_shortcut_creator()
