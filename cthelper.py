import ctypes as ct

# convert a ctypes character pointer to a string
def c_char_p2str(cp):
    return ct.cast(cp, ct.c_char_p).value

# convert a ctypes float pointer to a float
# N.B. one can inadvertantly reference off the end of the resulting list!!\
# it is best to trim based on the known length of the returned array!
def c_float_p2float(cp, numels):
    '''
    cp --> is a pointer to the start of a float array
    numel --> is the number of elements of the array to return
    ''' 
    return ct.cast(cp, ct.POINTER(ct.c_float))[0:numels]
    
# convert a ctypes double pointer to a float
# N.B. one can inadvertantly reference off the end of the resulting list!!\
# it is best to trim based on the known length of the returned array!
def c_double_p2float(cp, numels):
    '''
    cp --> is a pointer to the start of a double array
    numel --> is the number of elements of the array to return
    ''' 
    return ct.cast(cp, ct.POINTER(ct.c_double))[0:numels]
    