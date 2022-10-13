#Function for reading the matrix from the given text file
def read_Mat(file_name):
    with open(str(file_name)) as f:
        file = list(map(lambda x:list(map(lambda y:int(y),x)),
                        list(map(lambda x: x[1:-1].split(","),
                        str("".join(str("["+f.read()+"]").replace("\n","]@["))).split("@"))))) 
    return file
