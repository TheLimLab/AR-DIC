import http.server
import socketserver
import os
os.getcwd()
#set folder to serve:
os.chdir("C:\\Users\\username\\directory\\AR_DIC\\Synthetic_topography\\Tangram_mini_demo")
os.getcwd()
PORT = 8000

Handler = http.server.SimpleHTTPRequestHandler

with socketserver.TCPServer(("", PORT), Handler) as httpd:
    print("serving at port", PORT)
    httpd.serve_forever()
