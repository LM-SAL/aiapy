#Module for the Heliophysics Events Knowledgebase
#contributors: Nicholas Padmanabhan,  Mark Cheung

class HER_Event:
    def __init__(self,eventType):
        import numpy as np
        self.eventType = eventType #You need to eventually make eventType = 'CH: CoronalHole' instead of 'CH'
        self.specFile = '/ssw/vobs/ontology/data/VOEvent_Spec.txt'
        self.Reference_Names = np.zeros(shape=20,dtype=str)
        self.Reference_Links = np.zeros(shape=20,dtype=str)
        self.Reference_Types = np.zeros(shape=20,dtype=str)
        self.Description = ''
        self.Citations = np.zeros(shape=20,dtype=str)
        
        #str buffer to store XML code
        self.__xmlbuffer__ = ''

        # Some code to populate the required and optional dictionaries based on VOEvent_Spec
        self.required = {'kb_archivid':'Reserved for KB archivist: KB entry identifier','FRM_Name':'Dummy Name'}
        self.optional = {}
        # Nicholas to add code here
        
    def exportEvent(self,filename=None,suffix=''):
        import datetime

        # The AuthorIVORN is set by the FRM_Name and the time the export is performed
        remove_chars = [' ',':',';','/','\\', '{{','}}' ,'[[',']]', '~', '=', '.', ',', '-', '+']
        for r in remove_chars:
            self.required['FRM_Name'] = self.required['FRM_Name'].replace(r,'')

        #Set kb_archivid
        d = datetime.datetime.utcnow().isoformat() 
        eventIdentifier = (self.eventType.split(':'))[0] + '_' + self.required['FRM_Name'] + '_' + d + suffix
        self.required['kb_archivid'] = 'ivo://helio-informatics.org/' + eventIdentifier
        print(self.required['kb_archivid'])

        # Set filename if not given
        if (filename == None):
            filename = eventIdentifier+'.xml'
        
        # Some more lines of code to populate XML
        self.__xmlbuffer__ = "some placeholder text"
        # Nicholas to add code here
    
        # Write from buffer to file
        u = open(filename,"w+")
        u.write(self.__xmlbuffer__)
        u.close()