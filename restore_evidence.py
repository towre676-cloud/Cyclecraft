import pathlib,json,hashlib,tarfile,io
r=pathlib.Path(__file__).resolve().parent
idx=json.loads((r/"PAYLOAD_INDEX.json").read_text())
b=b"".join((r/x).read_bytes() for x in idx["parts"])
assert hashlib.sha256(b).hexdigest()==idx["sha256"]
with tarfile.open(fileobj=io.BytesIO(b),mode="r:gz") as t:
 for row in idx["files"]:
  data=t.extractfile(row["sha256"]).read()
  assert hashlib.sha256(data).hexdigest()==row["sha256"]
  p=r/row["path"];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)

import pathlib,json,hashlib
r=pathlib.Path(__file__).resolve().parent
for row in json.loads((r/"CHUNKED_FILES.json").read_text()):
 b=b"".join((r/x).read_bytes() for x in row["parts"])
 assert hashlib.sha256(b).hexdigest()==row["sha256"]
 p=r/row["path"];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)

import pathlib,json,gzip,hashlib
r=pathlib.Path(__file__).resolve().parent
for row in json.loads((r/"IMPORT_MANIFEST.json").read_text())["files"]:
 if row.get("storage")=="gzip":
  b=gzip.decompress((r/row["stored_path"]).read_bytes())
  assert hashlib.sha256(b).hexdigest()==row["sha256"], row["path"]
  p=r/row["path"];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
print("Restored and hash-verified all Git-hosted compressed evidence. See manifest for externally retained proof streams.")
