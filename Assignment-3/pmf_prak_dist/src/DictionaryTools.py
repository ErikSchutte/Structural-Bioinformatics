def find_key(dic, val):
    """return the key of dictionary dic given the value"""
    return [k for k, v in dic.iteritems() if v == val][0]

def find_value(dic, key):
    """return the value of dictionary dic given the key"""
    return dic[key]
