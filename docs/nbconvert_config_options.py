"""nbconvert config file"""

c = get_config() # noqa: F821
c.NbConvertApp.notebooks = ["../notebooks/**/*.ipynb"]
c.NbConvertApp.export_format = "rst"
c.NbConvertApp.recursive_glob = True
c.NbConvertApp.output_files_dir = "."
c.FilesWriter.build_directory = "."