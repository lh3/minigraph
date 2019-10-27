#!/usr/bin/env k8

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

var c, min_rm_div = 0.2, min_rm_sc = 300, micro_cap = 6, min_feat_len = 30;
while ((c = getopt(arguments, "")) != null) {
}

if (arguments.length - getopt.ind < 2) {
	print("Usage: anno-bb.js <bb.bed> <rm.out> [etrf.out] [sdust.out] [self.paf] [gap.bed]");
	exit(1);
}

var file, buf = new Bytes();

var bb = {}, bba = [];

file = new File(arguments[getopt.ind]);
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	var key = t[0] + "_" + t[1] + "_" + t[2];
	bb[key] = [parseInt(t[5]), {}, t[4], t[7], t[6], t[3]];
	bba.push(key);
}
file.close();

function process_rm_line(bb, lines)
{
	var h = {};
	if (lines.length == 0) return;
	var key = lines[0][4], h = bb[key][1];
	for (var i = 0; i < lines.length; ++i) {
		var t = lines[i];
		var st = parseInt(t[5]) - 1, en = parseInt(t[6]);
		if (h[t[10]] == null) h[t[10]] = [];
		h[t[10]].push([st, en]);
	}
}

file = new File(arguments[getopt.ind+1]);
var lines = [];
while (file.readline(buf) >= 0) {
	var line = buf.toString();
	var l2 = line.replace(/^\s+/, "");
	var t = l2.split(/\s+/);
	if (t.length < 15) continue;
	if (t[10] == 'Simple_repeat' || t[10] == 'Low_complexity') t[10] = 'LCR';
	if (t[10] != 'LCR') {
	//	if (parseInt(t[0]) < min_rm_sc) continue;
	//	if (parseInt(t[1])/100 > min_rm_div) continue;
	}
	if (lines.length > 0 && lines[0][4] != t[4]) {
		process_rm_line(bb, lines);
		lines = [];
	}
	lines.push(t);
}
if (lines.length > 0) process_rm_line(bb, lines);
file.close();

if (getopt.ind + 2 < arguments.length) {
	file = new File(arguments[getopt.ind+2]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var l = parseInt(t[4]);
		if (l == 1) continue;
		var anno = l <= micro_cap? 'micro' : 'mini';
		if (bb[t[0]][1][anno] == null)
			bb[t[0]][1][anno] = [];
		var st = parseInt(t[1]), en = parseInt(t[2]);
		bb[t[0]][1][anno].push([st, en]);
		if (bb[t[0]][1]['LCR'] == null)
			bb[t[0]][1]['LCR'] = [];
		bb[t[0]][1]['LCR'].push([st, en]);
	}
	file.close();
}

if (getopt.ind + 3 < arguments.length) {
	file = new File(arguments[getopt.ind+3]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var anno = 'LCR';
		if (bb[t[0]][1][anno] == null)
			bb[t[0]][1][anno] = [];
		bb[t[0]][1][anno].push([parseInt(t[1]), parseInt(t[2])]);
	}
	file.close();
}

if (getopt.ind + 4 < arguments.length) {
	file = new File(arguments[getopt.ind+4]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var anno = 'self';
		if (bb[t[0]][1][anno] == null)
			bb[t[0]][1][anno] = [];
		bb[t[0]][1][anno].push([parseInt(t[2]), parseInt(t[3])]);
	}
	file.close();
}

if (getopt.ind + 5 < arguments.length) {
	file = new File(arguments[getopt.ind+5]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var anno = 'gap';
		if (bb[t[0]][1][anno] == null)
			bb[t[0]][1][anno] = [];
		bb[t[0]][1][anno].push([parseInt(t[1]), parseInt(t[2])]);
	}
	file.close();
}

for (var i = 0; i < bba.length; ++i) {
	var m, key = bba[i], h = bb[key][1], len = bb[key][0];
	if ((m = /^(\S+)_(\d+)_(\d+)/.exec(key)) == null)
		throw("Bug!");
	var x = {}, t = [m[1], m[2], m[3], bb[key][0], bb[key][2], bb[key][3], bb[key][4], bb[key][5]];
	for (var c in h) {
		var s, st = 0, en = 0, cov = 0;
		s = h[c].sort(function(a, b) { return a[0] - b[0]; });
		for (var j = 0; j < s.length; ++j) {
			if (s[j][0] > en) {
				cov += en - st;
				st = s[j][0], en = s[j][1];
			} else en = en > s[j][1]? en : s[j][1];
		}
		cov += en - st;
		if (cov >= min_feat_len)
			x[c] = cov;
	}
	var type = "none";
	var max = 0, max2 = 0, max_c2 = null, max_c = null, sum = 0, sum_misc = 0;
	var lcr = x['LCR'] == null? 0 : x['LCR'];
	var self_len = x['self'] == null? 0 : x['self'];
	for (var c in x) {
		if (c == 'LCR' || c == 'self') continue;
		sum += x[c];
		if (c != 'mini' && c != 'micro') sum_misc += x[c];
		if (max < x[c]) max2 = max, max_c2 = max_c, max = x[c], max_c = c;
		else if (max2 < x[c]) max2 = x[c], max_c2 = c;
	}
	if (max >= len * 0.7) {
		type = max_c;
	} else if (lcr >= len * 0.7) {
		type = 'lcr';
		if (max_c == 'mini' || max_c == 'micro') {
			var y = x['mini'] == null? 0 : x['mini'];
			y += x['micro'] == null? 0 : x['micro'];
			if (max >= y * 0.7) type = max_c;
		}
	} else if ((max_c == 'mini' || max_c == 'micro') && max2 < max * 0.1) {
		type = max_c;
	} else if (sum_misc + lcr >= len * 0.7) {
		type = 'mixed';
	} else if (sum + lcr > len * 0.05) {
		type = 'partial';
	} else if (self_len >= len * 0.5) {
		type = 'self';
	}
	t.push(type);
	for (var c in x)
		t.push(c + ':' + x[c]);
	print(t.join("\t"));
}

buf.destroy();
