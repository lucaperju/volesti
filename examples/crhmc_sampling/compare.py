# written with chatgpt
import sys
import hashlib

def hash_file(filename):
    """This function returns the SHA-256 hash of the file passed into it"""
    h = hashlib.sha256()
    
    # Open the file in binary mode and update the hash
    with open(filename, 'rb') as file:
        while chunk := file.read(8192):
            h.update(chunk)
    
    return h.hexdigest()

def compare_files(file1, file2):
    """Compares two files by their hash values"""
    hash1 = hash_file(file1)
    hash2 = hash_file(file2)
    
    if hash1 == hash2:
        print(f"The files '{file1}' and '{file2}' are identical.")
    else:
        print(f"The files '{file1}' and '{file2}' are different.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    compare_files(file1, file2)
