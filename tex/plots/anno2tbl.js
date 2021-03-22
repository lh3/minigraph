#!/usr/bin/env k8

var buf = new Bytes();
var file = arguments.length == 0? new File() : new File(arguments[0]);

var h = {};
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	for (var i = 1; i <= 7; ++i) t[i] = parseInt(t[i]);
	if (t[5]) continue;
	if (t[11] == "gap") continue;
	if (/chrUn|_random/.test(t[0])) continue;
	var na = t[4] < 4? t[4] : 4;
	var type = null;
	if (t[11] == "mini") type = "11_VNTR";
	else if (t[11] == "micro") type = "12_STR";
	else if (t[11] == "micro" || t[11] == "lcr") type = "13_Other-LCR";
	else if (t[11] == "LINE/L1") type = "02_L1";
	else if (t[11] == "SINE/Alu") type = "01_Alu";
	else if (t[11] == "Retroposon/SVA") type = "03_SVA";
	else if (t[11] == "LTR/ERV") type = "04_ERV";
	else if (t[11] == "inter" || /^(DNA|LINE|SINE|LTR)/.test(t[11])) type = "05_Mixed-MEI";
	else if (/^Satellite/.test(t[11]) || t[11] == "alpha" || t[11] == "hsat2/3") type = "10_Satellite";
	else if (t[11] == "self") type = "31_Non-rep-dup";
	else if (t[11] == "none") type = "30_Non-rep-uniq";
	else if (t[11] == "mixed") type = "20_Mixed-repeat";
	else type = "21_Partial-repeat";
	var key = type;
	if (h[key] == null) h[key] = [0, null, 0, 0, 0, 0, 0, 0];
	++h[key][na];
	h[key][na+3] += t[7];
}

file.close();
buf.destroy();

for (var key in h) {
	var label = key.replace(/^[0-9]+_/, "");
	print(key, label, h[key][2], h[key][3], h[key][4], h[key][5], h[key][6], h[key][7]);
}
