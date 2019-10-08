import json
with open("CELL_LINES.json", "r") as f:
	d = json.loads(f.readline())
keys=sorted(d.keys())
n=0
for key in keys:
	if n ==2071:
		print(key)
		print(n)
	if key.find("/")!=-1:
		print(key)
		print(n)
	n+=1

