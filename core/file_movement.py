class FileMovement:
    def __init__(self,file_server):
        self.file_server = file_server
    def upload_file(self,configs):
        self.file_server.upload(configs)