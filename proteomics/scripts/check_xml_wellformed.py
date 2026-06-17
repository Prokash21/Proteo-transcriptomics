import sys
import xml.etree.ElementTree as ET
from pathlib import Path

p = Path(sys.argv[1]) if len(sys.argv)>1 else Path('PXD001676/Non-tumor-2_xml.xml')
print('Checking', p)
try:
    for i,(event,elem) in enumerate(ET.iterparse(p, events=['end'])):
        if i%500000==0:
            print('parsed events:', i)
    print('Parsed OK')
except ET.ParseError as e:
    print('PARSE ERROR:', e)
    sys.exit(2)
except Exception as e:
    print('ERROR:', e)
    sys.exit(3)
