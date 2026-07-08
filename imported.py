from imports import * #Imports

def check_imported_modules(modules_to_check): #Function to verify which modules are installed
    modules=sys.modules.keys() #Loads currently imported modules
    uninstalled_module=False #Parameter in case of missing modules
    counter=0 #Index to iterate along list of requested modules
    last_module_index=modules_to_check.index(modules_to_check[-1]) #Locates maximum necessary iterating index

    while uninstalled_module==False and counter<=last_module_index: #Loop that ensures all modules are checked and flagged as uninstalled if necessary
        for module in modules_to_check:
            if module not in modules: #Check for module installation
                print(f'{module} not found.') #Outputs name of missing module
                uninstalled_module=True #Flags module as missing
            counter+=1 #Iterates to next module
                