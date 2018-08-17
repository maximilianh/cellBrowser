# encode: utf-8

import cherrypy
WEB_ROOT = "/usr/local/var/www/"

class CServer( object ) :
    @cherrypy.expose
    def do_contact(self, **params):
        pass

cherrypy.server.socket_port = 8888
# INADDR_ANY: listen on all interfaces
#cherrypy.server.socket_host = '0.0.0.0'
cherrypy.server.socket_host = '127.0.0.1'

class Root(object): pass
conf = { '/':
  { 'tools.staticdir.on' : True,
    'tools.staticdir.dir' : WEB_ROOT,
    'tools.staticdir.index' : 'index.html' } }
#cherrypy.quickstart( CServer(), config = conf )
cherrypy.tree.mount(Root(), '/', config=conf)

cherrypy.quickstart()
