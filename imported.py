from imports import *

def check_imported_modules(modules_to_check):
    modules=sys.modules.keys()
    uninstalled_module=False
    counter = 0
    last_module_index=modules_to_check.index(modules_to_check[-1])

    while uninstalled_module == False and counter<=last_module_index:
        for module in modules_to_check:
            if module not in modules:
                print(f'{module} not found.')
                uninstalled_module=True
            counter+=1
                