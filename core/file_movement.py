class FileMovement:
    def __init__(self,file_server):
        self.file_server = file_server
    def upload_file(self,configs,file):
        self.file_server.upload(configs,file)