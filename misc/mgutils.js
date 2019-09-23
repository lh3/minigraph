#!/usr/bin/env k8

/*******************************
 * Command line option parsing *
 *******************************/

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

function it_index(a) {
	if (a.length == 0) return -1;
	a.sort(function(x, y) { return x[0] - y[0] });
	var last, last_i;
	for (var i = 0; i < a.length; i += 2) last = a[i][2] = a[i][1], last_i = i;
	for (var k = 1; 1<<k <= a.length; ++k) {
		var i0 = (1<<k) - 1, step = 1<<(k+1);
		for (var i = i0; i < a.length; i += step) {
			var x = 1<<(k-1);
			a[i][2] = a[i][1];
			if (a[i][2] < a[i-x][2]) a[i][2] = a[i-x][2];
			var e = i + x < a.length? a[i+x][2] : last;
			if (a[i][2] < e) a[i][2] = e;
		}
		last_i = last_i>>k&1? last_i - (1<<(k-1)) : last_i + (1<<(k-1));
		if (last_i < a.length) last = last > a[last_i][2]? last : a[last_i][2];
	}
	return k - 1;
}

function it_overlap(a, st, en) {
	if (a == null) return [];
	var h, stack = [], b = [];
	for (h = 0; 1<<h <= a.length; ++h);
	--h;
	stack.push([(1<<h) - 1, h, 0]);
	while (stack.length) {
		var t = stack.pop();
		var x = t[0], h = t[1], w = t[2];
		if (h <= 2) {
			var i0 = x >> h << h, i1 = i0 + (1<<(h+1)) - 1;
			if (i1 >= a.length) i1 = a.length;
			for (var i = i0; i < i1; ++i)
				if (a[i][0] < en && st < a[i][1])
					b.push(a[i]);
		} else if (w == 0) { // if left child not processed
			stack.push([x, h, 1]);
			var y = x - (1<<(h-1));
			if (y >= a.length || a[y][2] > st)
				stack.push([y, h - 1, 0]);
		} else if (x < a.length && a[x][0] < en) {
			if (st < a[x][1]) b.push(a[x]);
			stack.push([x + (1<<(h-1)), h - 1, 0]);
		}
	}
	return b;
}

function it_contained(a, st, en) {
	if (a == null) return false;
	var b = it_overlap(a, st, en);
	var c = false;
	for (var i = 0; i < b.length; ++i) {
		if (b[i][0] <= st && en <= b[i][1])
			c = true;
	}
	return c;
}

/****************************
 ***** mgutils commands *****
 ****************************/

function mg_cmd_renamefa(args)
{
	var c, sep = '#';
	while ((c = getopt(args, "d:")) != null)
		if (c == 'd') sep = getopt.arg;
	if (args.length - getopt.ind < 2) {
		print("Usage: mgutils.js renamefa [-d delimitor] <prefix> <in.fa>");
		return;
	}
	var prefix = args[getopt.ind];
	var file = new File(args[getopt.ind+1]);
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		if (buf[0] != 62) {
			print(buf);
		} else {
			var m, s = buf.toString();
			if ((m = /^>(.*)/.exec(s)) != null) {
				var name = m[1].replace(/^\S+#/, "");
				print(">" + prefix + sep + name);
			} else throw Error("Wrong FASTA format!");
		}
	}
	file.close();
	buf.destroy();
}

function mg_cmd_subgaf(args)
{
	if (args.length < 2) {
		print("Usage: mgutils.js subgaf <in.gaf> <reg>");
		exit(1);
	}

	var m, ctg, st, en;
	if ((m = /^(\S+):(\S+)-(\S+)/.exec(args[1])) != null)
		ctg = m[1], st = parseInt(m[2]), en = parseInt(m[3]);

	var buf = new Bytes();
	var file = new File(args[0]);
	var re = /([><])([^\s><]+):(\d+)-(\d+)/g;

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var l = parseInt(t[6]), s = parseInt(t[7]), e = parseInt(t[8]);
		var regs = [];
		if (t[5][0] == '>' || t[5][0] == '<') {
			var m, x = 0;
			//print(buf);
			while ((m = re.exec(t[5])) != null) {
				var a = parseInt(m[3]), b = parseInt(m[4]), c = b - a;
				if (x == 0) {
					if (b - a <= s) throw Error("Inconsistent!");
					a += s;
				}
				if (x + c == l) b -= l - e;
				//print(m[2], a, b);
				regs.push([m[2], a, b]);
				x += c;
			}
		} else {
			regs.push([t[5], s, e]);
		}
		var hit = false;
		for (var i = 0; i < regs.length; ++i) {
			if (regs[i][0] == ctg && regs[i][2] > st && en > regs[i][1])
				hit = true;
		}
		if (hit) print(buf);
	}

	file.close();
	buf.destroy();
}

function mg_cmd_sveval(args)
{
	var c, flank = 100, min_var_len = 100, min_test_len = 50, min_sc = 20.0, non_chr = false, out_err = false;
	while ((c = getopt(args, "f:v:t:s:ae")) != null) {
		if (c == 'f') flank = parseInt(getopt.arg);
		else if (c == 'v') min_var_len = parseInt(getopt.arg);
		else if (c == 't') min_test_len = parseInt(getopt.arg);
		else if (c == 's') min_sc = parseFloat(getopt.arg);
		else if (c == 'a') non_chr = true;
		else if (c == 'e') out_err = true;
	}
	if (args.length - getopt.ind < 3) {
		print("Usage: mgutils.js sveval <true.vcf> <true.bed> <call.txt>");
		print("Options:");
		print("  -f INT      length of flanking regions [" + flank + "]");
		print("  -v INT      min INDEL length [" + min_var_len + "]");
		print("  -t INT      min true INDEL length [" + min_test_len + "]");
		print("  -s INT      min called score [" + min_sc + "]");
		print("  -e          print errors");
		exit(1);
	}

	var file, buf = new Bytes();

	// parse true.bed
	warn("Reading confident regions...");
	var bed = {}
	file = new File(args[getopt.ind + 1]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t.length < 3) continue;
		if (!non_chr && /^(chr)?[XY]$/.test(t[0])) continue;
		if (bed[t[0]] == null) bed[t[0]] = [];
		bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	file.close();
	for (var ctg in bed) it_index(bed[ctg]);

	// parse true.vcf
	warn("Reading baseline variants...");
	var vcf = {}, n_vcf = 0;
	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		if (t.length < 10) continue;
		if (t[6] != '.' && t[6] != 'PASS') continue;
		if (bed[t[0]] == null) continue;
		var ref = t[3];
		var st = parseInt(t[1]) - 1;
		var en = st + ref.length;
		var max_diff = 0;
		var al = t[4].split(",");
		al.unshift(ref);
		for (var i = 1; i < al.length; ++i) {
			var l = al[i].length - ref.length;
			if (l < 0) l = -l;
			if (max_diff < l) max_diff = l;
		}
		if (max_diff < min_test_len) continue;
		var s = t[9].split(':');
		if (s.length == 0) continue;
		var gt = s[0].split(/[|\/]/);
		if (gt == 0) continue;
		var max_ev = 0;
		max_diff = 0;
		for (var i = 0; i < gt.length; ++i) {
			var x = parseInt(gt[i]);
			var l = al[x].length - ref.length;
			var x = l > 0? l : -l;
			if (max_diff < x) max_diff = x, max_ev = l;
		}
		if (max_diff < min_test_len) continue;
		if (vcf[t[0]] == null) vcf[t[0]] = [];
		vcf[t[0]].push([st, en, -1, max_diff, max_ev]);
	}
	file.close();
	for (var ctg in vcf) it_index(vcf[ctg]);

	// parse rst.txt
	warn("Reading gt results...");
	var rst = {};
	file = new File(args[getopt.ind + 2]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (parseFloat(t[3]) < min_sc) continue;
		if (bed[t[0]] == null) continue;
		if (rst[t[0]] == null) rst[t[0]] = [];
		var st = parseInt(t[1]), en = parseInt(t[2]);
		rst[t[0]].push([st, en]);
	}
	file.close();
	for (var ctg in rst) it_index(rst[ctg]);

	// sensitivity
	var n_vcf = [0, 0, 0], fn = [0, 0, 0];
	for (var ctg in vcf) {
		for (var i = 0; i < vcf[ctg].length; ++i) {
			var v = vcf[ctg][i];
			if (!it_contained(bed[ctg], v[0], v[1])) continue;
			if (v[3] < min_var_len) continue;
			var sub = v[4] < 0? 1 : 2;
			++n_vcf[0], ++n_vcf[sub];
			var st = v[0] - flank, en = v[1] + flank;
			if (st < 0) st = 0;
			var b = it_overlap(rst[ctg], st, en);
			if (b.length == 0) {
				if (out_err) print("FN", ctg, v[0], v[1], v[4]);
				++fn[0], ++fn[sub];
			}
		}
	}

	// specificity
	var n_rst = 0, fp = 0;
	for (var ctg in rst) {
		for (var i = 0; i < rst[ctg].length; ++i) {
			var v = rst[ctg][i];
			if (!it_contained(bed[ctg], v[0], v[1])) continue;
			++n_rst;
			var st = v[0] - flank, en = v[1] + flank;
			if (st < 0) st = 0;
			var b = it_overlap(vcf[ctg], st, en);
			if (b.length == 0) {
				if (out_err) print("FP", ctg, v[0], v[1]);
				++fp;
			}
		}
	}

	print("NA", fn[0], n_vcf[0], (fn[0]/n_vcf[0]).toFixed(4));
	print("ND", fn[1], n_vcf[1], (fn[1]/n_vcf[1]).toFixed(4));
	print("NI", fn[2], n_vcf[2], (fn[2]/n_vcf[2]).toFixed(4));
	print("PA", fp, n_rst, (fp/n_rst).toFixed(4));
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: mgutils.js <command> [arguments]");
		print("Commands:");
		print("  renamefa     add a prefix to sequence names in FASTA");
		print("  subgaf       extract GAF overlapping with a region");
		print("  sveval       evaluate SV accuracy");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'renamefa') mg_cmd_renamefa(args);
	else if (cmd == 'subgaf') mg_cmd_subgaf(args);
	else if (cmd == 'sveval') mg_cmd_sveval(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
