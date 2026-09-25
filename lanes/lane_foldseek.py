"""Foldseek API lane: structural homology search of the 8A8C docking-control
receptor against the PDB (validates the AlphaFold/receptor structure set)."""
import json, time, urllib.request, urllib.parse

def req(url, data=None, files=None):
    r = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(r, timeout=60) as fh:
        return json.load(fh)

def main():
    pdb = open("data/raw/structures/8A8C.pdb").read()
    boundary = "----mega27boundary"
    parts = []
    for name, val in (("q", pdb), ("database[]", "pdb100"), ("mode", "3diaa")):
        parts.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\""
                     + ("; filename=\"q.pdb\"\r\nContent-Type: chemical/x-pdb" if name == "q" else "")
                     + f"\r\n\r\n{val}\r\n")
    body = ("".join(parts) + f"--{boundary}--\r\n").encode()
    r = urllib.request.Request("https://search.foldseek.com/api/ticket", data=body,
                               headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(r, timeout=60) as fh:
        ticket = json.load(fh)
    tid = ticket["id"]
    status = ""
    for _ in range(40):
        time.sleep(6)
        st = req(f"https://search.foldseek.com/api/ticket/{tid}")
        status = st.get("status", "")
        if status in ("COMPLETE", "ERROR"):
            break
    out = {"ticket": tid, "status": status, "hits": []}
    if status == "COMPLETE":
        res = req(f"https://search.foldseek.com/api/result/{tid}/0")
        for e in res.get("results", [{}])[0].get("alignments", [])[:15]:
            out["hits"].append({"target": e.get("target"), "seqid": e.get("seqId"),
                                "alnlen": e.get("alnLength"), "evalue": e.get("eval")})
    json.dump(out, open("results/foldseek_8a8c_pdb_search.json", "w"), indent=2)
    print(status, len(out["hits"]))

if __name__ == "__main__":
    main()
