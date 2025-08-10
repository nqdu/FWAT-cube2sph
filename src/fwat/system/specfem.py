def change_parfile(filename,**kargs):
    with open(filename,"r") as fio:
        lines = fio.readlines()

    # loop each value in dictionary to set in lines
    for key,val in kargs.items():
        nline = len(lines)
        find_it = False
        for i in range(nline):
            if len(lines[i]) <= 1 or lines[i][0] =='#':
                continue

            keystr = key.upper() + " "
            n = len(keystr)
            if keystr == lines[i][:n]:
                # find it
                find_it = True
                lines[i] = f"{key}         =  {val}\n"
                break 
    
        # check if value in it
        if not find_it:
            print(f"{key} donnot exist in Par_file")
            exit(1)
    if filename is None:
        return '\n'.join(lines)
    else:
        with open(filename,"w") as fio:
            fio.writelines(lines)

def get_param(filename:str,key:str):
    with open(filename,"r") as fio:
        lines = fio.readlines()

    nline = len(lines)
    for i in range(nline):
        if len(lines[i]) <= 1 or lines[i][0] =='#':
            continue
        
        var = ""
        keystr = key.upper() + " "
        n = len(keystr)
        if keystr == lines[i][:n]:
            # find it
            find_it = True

            # get it 
            val = lines[i].split('=')[-1].split()[0]

            break 

    # check if value in it
    if not find_it:
        print(f"{key} donnot exist in Par_file")
        exit(1)

    return val