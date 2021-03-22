#!/usr/bin/env k8

/*****************************
 ***** Library functions *****
 *****************************/

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

/***************
 * BED overlap *
 ***************/

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
					b.push(i);
		} else if (w == 0) { // if left child not processed
			stack.push([x, h, 1]);
			var y = x - (1<<(h-1));
			if (y >= a.length || a[y][2] > st)
				stack.push([y, h - 1, 0]);
		} else if (x < a.length && a[x][0] < en) {
			if (st < a[x][1]) b.push(x);
			stack.push([x + (1<<(h-1)), h - 1, 0]);
		}
	}
	return b;
}

/******************************
 ***** Command-line tools *****
 ******************************/

function bed_sum(args)
{
	var buf = new Bytes();
	var file = args.length == 0 || args[0] == '-'? new File() : new File(args[0]);
	var s = 0;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 3);
		if (t.length < 3) continue;
		s += parseInt(t[2])  - parseInt(t[1]);
	}
	file.close();
	buf.destroy();
	print(s);
	return 0;
}

function bed_sum2nd(args)
{
	var buf = new Bytes();
	var file = args.length == 0 || args[0] == '-'? new File() : new File(args[0]);
	var s = 0;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 2);
		s += parseInt(t[1]);
	}
	file.close();
	buf.destroy();
	print(s);
	return 0;
}

function bed_merge(args)
{
	var buf = new Bytes();
	var file = args.length > 0? new File(args[0]) : new File();
	var ctg = null, st, en;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 3);
		var s = parseInt(t[1]);
		var e = parseInt(t[2]);
		if (ctg != t[0] || s > en) { // no overlap
			if (ctg != null) print(ctg, st, en);
			ctg = t[0], st = s, en = e;
		} else if (s < st) throw Error("ERROR: input is not sorted by coordinate");
		else en = en > e? en : e;
	}
	if (ctg != null) print(ctg, st, en);
	file.close();
	buf.destroy();
	return 0;
}

function bed_sum1(args)
{
	var buf = new Bytes();
	var file = args.length == 0 || args[0] == '-'? new File() : new File(args[0]);
	var ctg = null, st = 0, en = 0, sum = 0;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 3);
		var s = parseInt(t[1]);
		var e = parseInt(t[2]);
		if (ctg != t[0] || s > en) { // no overlap
			sum += en - st;
			if (ctg != null && ctg != t[0]) {
				print(ctg, sum);
				sum = 0;
			}
			ctg = t[0], st = s, en = e;
		} else if (s < st) throw Error("ERROR: input is not sorted by coordinate");
		else en = en > e? en : e;
	}
	if (ctg != null) {
		sum += en - st;
		print(ctg, sum);
	}
	file.close();
	buf.destroy();
	return 0;
}

function bed_gdist(args)
{
	if (args.length == 0) {
		print("Usage: bedutils.js gdist <3-col-gmap.txt> <reg.bed>");
		exit(1);
	}
	var file, buf = new Bytes();

	var gmap = {};
	file = new File(args[0]);
	var last_pos = 0, last_ctg = null, last_v = 0.0;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var pos = parseInt(t[1]);
		var v = parseFloat(t[2]);
		if (last_ctg != t[0] && last_ctg != null) {
			gmap[last_ctg].push([last_pos, 0x7fffffff, -1, last_v]);
			last_pos = 0, last_v = 0.0;
		}
		if (gmap[t[0]] == null) gmap[t[0]] = [];
		if (last_pos == pos) throw Error("Zero-length interval");
		gmap[t[0]].push([last_pos, pos, -1, last_v]);
		last_pos = pos, last_ctg = t[0], last_v = v;
	}
	if (last_ctg != null)
		gmap[last_ctg].push([last_pos, 0x7fffffff, -1, last_v]);
	file.close();

	for (var ctg in gmap) it_index(gmap[ctg]);

	file = args.length >= 2? new File(args[1]) : new File();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var st = parseInt(t[1]), en = parseInt(t[2]);
		var v, g = gmap[t[0]];
		if (g == null) v = -1;
		else if (st == en) v = 0;
		else {
			var as = it_overlap(g, st, st + 1);
			var ae = it_overlap(g, en - 1, en);
			if (as.length != 1 || ae.length != 1)
				throw Error("Bug!");
			var is = as[0], ie = ae[0];
			var xs = g[is][3] + (is == g.length - 1? 0 : (g[is+1][3] - g[is][3]) / (g[is][1] - g[is][0]) * (st - g[is][0]));
			var xe = g[ie][3] + (ie == g.length - 1? 0 : (g[ie+1][3] - g[ie][3]) / (g[ie][1] - g[ie][0]) * (en - g[ie][0]));
			v = 1e6 * (xe - xs) / (en - st);
		}
		v = v <= 0? v : v.toFixed(15);
		print(t[0], t[1], t[2], v);
	}
	file.close();
	buf.destroy();
}

function bed_window(args)
{
	var c, win_size = 1000000, skip = 500000, cnt_only = false, fn_len = null;
	while ((c = getopt(args, "w:s:cl:")) != null) {
		if (c == 'w') win_size = parseInt(getopt.arg);
		else if (c == 's') skip = parseInt(getopt.arg);
		else if (c == 'c') cnt_only = true;
		else if (c == 'l') fn_len = getopt.arg;
	}

	var lens = {}, file, buf = new Bytes();
	if (fn_len) {
		file = new File(fn_len);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (t.length < 2) continue;
			lens[t[0]] = parseInt(t[1]);
		}
		file.close();
	}
	file = getopt.ind < args.length? new File(args[getopt.ind]) : new File();
	var bed = {}, ctgs = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (bed[t[0]] == null) { bed[t[0]] = []; ctgs.push(t[0]); }
		bed[t[0]].push([parseInt(t[1]), parseInt(t[2]), -1]);
	}
	file.close();
	buf.destroy();

	for (var ct = 0; ct < ctgs.length; ++ct) {
		var ctg = ctgs[ct];
		it_index(bed[ctg]);
		var a = bed[ctg];
		var max = 0;
		for (var i = 0; i < a.length; ++i)
			max = max > a[i][1]? max : a[i][1];
		if (lens[ctg] > 0 && max < lens[ctg]) max = lens[ctg];
		for (var x = 0; x < max; x += skip) {
			var st = x - (win_size>>1), en = x + (win_size>>1);
			if (st < 0) st = 0;
			if (en > max) en = max;
			var sum = 0, b = it_overlap(a, st, en);
			if (cnt_only) {
				sum = b.length;
			} else {
				for (var i = 0; i < b.length; ++i) {
					var c = a[b[i]];
					var s = st > c[0]? st : c[0];
					var e = en < c[1]? en : c[1];
					sum += e - s;
				}
			}
			print(ctg, x, sum/(en-st)*1e6);
		}
	}
}

function bed_cov(args)
{
	if (args.length < 2) {
		warn("Usage: bedutils.js cov <loaded.bed> <streamed.bed>");
		exit(1);
	}
	var file, buf = new Bytes();

	file = new File(args[0]);
	var bed = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 3);
		if (bed[t[0]] == null) bed[t[0]] = [];
		bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	for (var ctg in bed) it_index(bed[ctg]);
	file.close();

	file = new File(args[1]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t", 3);
		if (bed[t[0]] == null) {
			print(t[0], t[1], t[2], 0, 0);
		} else {
			var st0 = parseInt(t[1]), en0 = parseInt(t[2]);
			var b = bed[t[0]];
			var a = it_overlap(b, st0, en0);
			var cov_st = 0, cov_en = 0, cov = 0;
			for (var i = 0; i < a.length; ++i) {
				var st1 = b[a[i]][0] > st0? b[a[i]][0] : st0;
				var en1 = b[a[i]][1] < en0? b[a[i]][1] : en0;
				if (st1 > cov_en) {
					cov += cov_en - cov_st;
					cov_st = st1, cov_en = en1;
				} else cov_en = cov_en > en1? cov_en : en1;
			}
			cov += cov_en - cov_st;
			print(t[0], t[1], t[2], a.length, cov);
		}
	}
	file.close();

	buf.destroy();
}

function main(args)
{
	if (args.length == 0) {
		print("Usage: bedutils.js <command> [arguments]");
		print("Commands:");
		print("  sum        sum of BED regions (deprecated by bedtk)");
		print("  sum1       sum of BED regions for each contig");
		print("  sum2nd     sum of the 2nd column");
		print("  merge      merge overlapping regions in *sorted* BED (deprecated)");
		print("  cov        breadth of coverage (deprecated by bedtk)");
		print("  gdist      genetic distance from 3-col genetic map");
		print("  window     window-based counting");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'sum') bed_sum(args);
	else if (cmd == 'sum2nd') bed_sum2nd(args);
	else if (cmd == 'sum1') bed_sum1(args);
	else if (cmd == 'merge') bed_merge(args);
	else if (cmd == 'cov') bed_cov(args);
	else if (cmd == 'gdist') bed_gdist(args);
	else if (cmd == 'window') bed_window(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
