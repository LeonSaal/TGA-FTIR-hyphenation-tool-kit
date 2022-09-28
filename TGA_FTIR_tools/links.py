
import logging
from os import PathLike

import requests

logger = logging.getLogger(__name__)


class LINKS():
    API = 'https://api.github.com/repos/LeonSaal/TGA-FTIR-hyphenation-tool-kit/contents'
    WIKI = 'https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/wiki'
    REPO = 'https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit#readme'

def download_supplementary(directory: str, filename: str, dst: str|PathLike):
    base_url = f'{LINKS.API}/{directory}'
    resp = requests.get(base_url)
    if resp.ok:
        contents = resp.json()
        for file in contents:
            if file['name']==filename:
                download_url = file['download_url']
                if (resp:=requests.get(download_url)).ok:
                    with open(dst, "wb") as file:
                        file.write(resp.content)
                        logger.info(f"Downloaded {filename!r} from repository and stored in {dst!r}.")
                break
    
    else:
        logger.error(f"Unable to download {filename!r} from repository.")
