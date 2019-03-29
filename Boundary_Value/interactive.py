import sys

def readArgs(args, params):
    """
    Function which translates command line arguments and acts on them.
    :param args: List of arguments to the command line.
    """
    i=0
    updates = {}
    while i < len(args):
        if args[i] in ["help", "Help", "-h", "H", "h", "?", "??"]:
            with open("README.md", "r") as readme:
                print(readme.read())
                exit()
        elif args[i] in ["-RS", "-rs"]: 
            try:
                updates["Seed"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised seed.")
                exit()
        elif args[i] in ["-x", "-X"]:
            try:
                updates["X Dimension"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -x.")
                exit()
        elif args[i] in ["-D", "-d"]:
            try:
                if args[i+1] == "J":
                    updates["Dynamics"] = "Jacobi"
                elif args[i+1] in ["G", "GS"]:
                    updates["Dynamics"] = "GS"
                i += 2
            except:
                print("Error with -d tag.")
        elif args[i] in ["-N", "-n"]:
            try:
                updates["tMax"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -N.")
                exit()
        elif args[i] in ["-y", "-Y"]:
            try:
                updates["Y Dimension"] = int(float(args[i+1]))
                i += 2
            except:
                print("Unrecognised value for -y.")
                exit()
        elif args[i] in ["-S", "-s"]:
            try:
                if args[i+1] in ["Y", "y"]:
                    updates["Animate"] = True
                    i += 2
                elif args[i+1] in ["N", "n"]:
                    updates["Animate"] = False
                    i += 2
                else:
                    print("-s should be followed by 'Y' or 'N'.")
                    exit()
            except:
                print("Error with -S tag.")
                exit()
        elif args[i] in ["-i", "-I"]:
            try:
                updates["phi0"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -I")
                exit()  
        elif args[i] in ["-w", "-W"]:
            try:
                updates["omega"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -W")
                exit()  
        elif args[i] in ["-t", "-T"]:
            try:
                updates["tolerance"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -T")
                exit() 
        elif args[i] in ["-a", "-A"]:
            try:
                updates["a"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -a")
                exit()   
        elif args[i] in ["-b", "-B"]:
            try:
                updates["b"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -b")
                exit()                
        elif args[i] in ["-K", "-k"]:
            try:
                updates["kappa"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -k")
                exit()
        elif args[i] in ["-m", "-M"]:
            try:
                updates["M"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -M")
                exit()
        elif args[i] in ["-u", "-U"]:
            try:
                updates["Rate"] = float(args[i+1])
                i += 2
            except:
                print("Unrecognised value for -U")
                exit()
        elif args[i] in ["-p", "-P"]:
            try:
                if args[i+1] in ["Y", "y"]:
                    updates["Measure"] = True
                    i += 2
                elif args[i+1] in ["N", "n"]:
                    updates["Measure"] = False
                    i += 2
                else:
                    print("-P should be followed by 'Y' or 'N'.")
                    exit()
            except:
                print("Error with -P tag.")
                exit()
        elif args[i] in ["-r", "-R"]:
            try:
                updates["RunLabel"] = args[i+1]
                i += 2
            except:
                print("Unrecognised value for -r.")
                exit()
        elif args[i] in ["-o", "-O"]:
            try:
                updates["outDir"] = args[i+1]
                i += 2
            except:
                print("Unrecognised value for -o.")
                exit()
        else:
            print("Key {} not recognised. Ignoring.".format(args[i]))
            i += 2
    params.update(updates)
        
        
                
            
            
