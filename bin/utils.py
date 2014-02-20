import os
from os import path
import httplib2
import urllib
import json

def check_dir(dirname):
    if not path.exists(dirname):
        os.mkdir(dirname)

# Ensembl REST API stuff
http = httplib2.Http(".cache")

server = "http://beta.rest.ensembl.org/"
def ens_get(ext, *args, **kwargs):
    if len(args):
        ext += '&'.join([ urllib.quote(a) for a in args]) 
        
    if len(kwargs):
        ext += '&' + urllib.urlencode(kwargs)

    print ext

    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    
    if not resp.status == 200:
        exc = IOError(resp)
        exc.errno = resp
        raise exc

    decoded = json.loads(content)
    return decoded


