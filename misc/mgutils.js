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

function mg_cmd_joinfa(args)
{
	var c, len_n = 20, min_len = 150, name = "decoy-cat";
	while ((c = getopt(args, "n:l:s:")) != null) {
		if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'n') len_n = parseInt(getopt.arg);
		else if (c == 's') name = getopt.arg;
	}
	if (args.length - getopt.ind < 1) {
		print("Usage: mgutils.js joinfa [options] <in.fa>");
		return;
	}
	var seq = new Bytes(), seq1 = new Bytes(), lineno = 0, nn = new Bytes();
	for (var i = 0; i < len_n; ++i) nn.set(78);
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		++lineno;
		if (buf[0] == 62) {
			if (seq1.length >= min_len) {
				if (seq.length > 0) seq.set(nn);
				seq.set(seq1);
			}
			seq1.length = 0;
		} else seq1.set(buf);
	}
	if (seq1.length >= min_len) {
		if (seq.length > 0) seq.set(nn);
		seq.set(seq1);
	}
	print(">" + name);
	print(seq);
	file.close();
	buf.destroy();
	seq.destroy();
	seq1.destroy();
}

function mg_cmd_anno(args)
{
	var c, min_rm_div = 0.2, min_rm_sc = 300, micro_cap = 6, min_feat_len = 30, min_centro_len = 200, mobile = false, max_mobile_div = 2.0;
	var fn_rmout = null, fn_etrf = null, fn_dust = null, fn_gap = null, fn_paf = null, fn_centro = null, fn_bb = null, fn_sd = null;
	while ((c = getopt(args, "e:p:g:d:r:c:l:b:s:m")) != null) {
		if (c == 'l') min_feat_len = parseInt(getopt.arg);
		else if (c == 'm') mobile = true;
		else if (c == 'e') fn_etrf = getopt.arg;
		else if (c == 'p') fn_paf = getopt.arg;
		else if (c == 'g') fn_gap = getopt.arg;
		else if (c == 'd') fn_dust = getopt.arg;
		else if (c == 'r') fn_rmout = getopt.arg;
		else if (c == 'c') fn_centro = getopt.arg;
		else if (c == 'b') fn_bb = getopt.arg;
		else if (c == 's') fn_sd = getopt.arg;
	}

	if (args.length - getopt.ind < 1) {
		print("Usage: anno.js [options] <in.bed>");
		print("Options:");
		print("  -l INT      min feature length [" + min_feat_len + "]");
		print("  -r FILE     RepeatMasker .out [null]");
		print("  -g FILE     seqtk gap output for stretches of Ns [null]");
		print("  -d FILE     minimap2/sdust output for LCRs [null]");
		print("  -e FILE     etrf output [null]");
		print("  -p FILE     PAF alignment against reference [null]");
		print("  -c FILE     dna-brnn centromere results [null]");
		print("  -b FILE     bubble file [null]");
		print("  -s FILE     segdup file (paste gfa2bed bedcov) [null]");
		print("  -m          annotate AluY and L1HS separately");
		exit(1);
	}

	var file, buf = new Bytes();

	var bb = {}, bba = [], seg = {};

	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t.length < 4) continue;
		var key = t[0] + "_" + t[1] + "_" + t[2];
		var len = parseInt(t[3]);
		if (len < parseInt(t[2]) - parseInt(t[1]))
			throw Error("ERROR: event length smaller than interval length");
		bb[key] = [len, {}];
		bba.push(key);
	}
	file.close();

	if (fn_bb) {
		if (fn_sd) {
			file = new File(fn_sd);
			while (file.readline(buf) >= 0) {
				var t = buf.toString().split("\t");
				seg[t[3]] = [parseInt(t[4]), parseInt(t[2]) - parseInt(t[1]), parseInt(t[6])];
			}
			file.close();
		}
		file = new File(fn_bb);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var key = t[0] + "_" + t[1] + "_" + t[2];
			if (key in bb) {
				bb[key].push(t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10]);
				var s = t[11].split(","), tot_len = 0, tot_sd = 0, ref_len = 0;
				for (var i = 1; i < s.length - 1; ++i) {
					if (seg[s[i]] == null) continue;
					tot_len += seg[s[i]][1], tot_sd += seg[s[i]][2];
					if (seg[s[i]][0] == 0)
						ref_len += seg[s[i]][1];
				}
				bb[key][7] = tot_len;
				bb[key][8] = tot_sd;
				bb[key][9] = ref_len;
			}
		}
		file.close();
	}

	if (fn_rmout) {
		var motif0 = "GGAAT", motif_hash = {}, motif_mut_hash = {};
		{ // dealing with possible (GGAAT)n rotations and mutations
			var comp_tbl = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' };
			var motif = [motif0], motif_alt = [];

			// reverse complement
			for (var i = 0; i < motif.length; ++i) {
				var x = motif[i], y = "";
				for (var j = x.length - 1; j >= 0; --j) {
					y += comp_tbl[x[j]];
				}
				motif_alt.push(y);
			}
			for (var i = 0; i < motif_alt.length; ++i)
				motif.push(motif_alt[i]);

			// rotate
			motif_alt = [];
			for (var i = 0; i < motif.length; ++i) {
				var x = motif[i];
				for (var j = 1; j < x.length; ++j)
					motif_alt.push(x.substr(j) + x.substr(0, j));
			}
			for (var i = 0; i < motif_alt.length; ++i)
				motif.push(motif_alt[i]);

			for (var i = 0; i < motif.length; ++i) motif_hash[motif[i]] = i;

			// mutate
			var bases = [ 'A', 'C', 'G', 'T' ];
			for (var x in motif_hash) {
				var y = x;
				for (var i = 0; i < x.length; ++i) {
					for (var j = 0; j < bases.length; ++j) {
						var a = x.split("");
						if (a[i] == bases[j]) continue;
						a[i] = bases[j];
						motif_mut_hash[a.join("")] = 1;
					}
				}
			}
		}

		function process_rm_line(bb, lines) {
			var h = {};
			if (lines.length == 0) return;
			var key = lines[0][4];
			if (bb[key] == null) throw Error("ERROR: missing key: " + key);
			var h = bb[key][1];
			for (var i = 0; i < lines.length; ++i) {
				var t = lines[i];
				var st = parseInt(t[5]) - 1, en = parseInt(t[6]);
				if (h[t[10]] == null) h[t[10]] = [];
				h[t[10]].push([st, en]);
			}
		}

		file = new File(fn_rmout);
		var lines = [];
		while (file.readline(buf) >= 0) {
			var line = buf.toString();
			var l2 = line.replace(/^\s+/, "");
			var m4, t = l2.split(/\s+/);
			if (t.length < 15) continue;
			if (t[9] == "ALR/Alpha") t[10] = "alpha";
			else if (t[9] == "HSATII") t[10] = "hsat2/3";
			else if (/^LTR\/ERV/.test(t[10])) t[10] = 'LTR/ERV';
			else if (/^LTR/.test(t[10])) t[10] = 'LTR/misc';
			else if (/^DNA/.test(t[10])) t[10] = 'DNA/misc';
			else if (/rRNA|scRNA|snRNA|srpRNA/.test(t[10])) t[10] = 'RNAmisc';
			else if (/^LINE/.test(t[10]) && t[10] != "LINE/L1") t[10] = 'LINE/misc';
			else if ((t[10] == "Simple_repeat" || t[10] == "Satellite") && ((m4 = /^\(([ACGT]+)\)n/.exec(t[9])) != null)) {
				if (motif_hash[m4[1]] != null) {
					t[10] = "hsat2/3";
				} else if (m4[1].length % motif0.length == 0) {
					var c = 0, c_mut = 0;
					for (var j = 0; j < m4[1].length; j += motif0.length) {
						var s = m4[1].substr(j, j + motif0.length);
						if (motif_hash[s] != null)
							++c;
						else if (motif_mut_hash[s] != null)
							++c_mut;
					}
					if (c > 0 && (c + c_mut) * motif0.length == m4[1].length)
						t[10] = "hsat2/3";
				}
			}

			if (mobile) {
				if (t[10] == "LINE/L1" && t[9] == "L1HS" && parseFloat(t[1]) < max_mobile_div) t[10] = "LINE/L1HS";
				if (t[10] == "SINE/Alu" && /^AluY/.test(t[9]) && parseFloat(t[1]) < max_mobile_div) t[10] = "SINE/AluY";
			}
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

		for (var i = 0; i < bba.length; ++i) {
			var h = bb[bba[i]][1], a = [];
			for (var key in h)
				if (/^(DNA|SINE|LINE|Retroposon|LTR)/.test(key))
					for (var j = 0; j < h[key].length; ++j)
						a.push(h[key][j]);
			if (a.length) h['_inter'] = a;
		}
	}

	if (fn_etrf) {
		file = new File(fn_etrf);
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

	if (fn_dust) {
		file = new File(fn_dust);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var anno = 'LCR';
			if (bb[t[0]][1][anno] == null)
				bb[t[0]][1][anno] = [];
			bb[t[0]][1][anno].push([parseInt(t[1]), parseInt(t[2])]);
		}
		file.close();
	}

	if (fn_paf) {
		file = new File(fn_paf);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var anno = 'self';
			if (bb[t[0]][1][anno] == null)
				bb[t[0]][1][anno] = [];
			bb[t[0]][1][anno].push([parseInt(t[2]), parseInt(t[3])]);
		}
		file.close();
	}

	if (fn_gap) {
		file = new File(fn_gap);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var anno = 'gap';
			if (bb[t[0]][1][anno] == null)
				bb[t[0]][1][anno] = [];
			bb[t[0]][1][anno].push([parseInt(t[1]), parseInt(t[2])]);
		}
		file.close();
	}

	if (fn_centro) { // this not used, as RepeatMasker can usually do a good job
		file = new File(fn_centro);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var anno = t[3] == '1'? 'hsat2/3' : 'alpha';
			if (bb[t[0]][1][anno] == null)
				bb[t[0]][1][anno] = [];
			var st = parseInt(t[1]), en = parseInt(t[2]);
			if (en - st >= min_centro_len)
				bb[t[0]][1][anno].push([st, en]);
		}
		file.close();
	}

	for (var i = 0; i < bba.length; ++i) {
		var m, key = bba[i], h = bb[key][1], len = bb[key][0];
		if ((m = /^(\S+)_(\d+)_(\d+)/.exec(key)) == null)
			throw("Bug!");
		var x = {}, t = [m[1], m[2], m[3]];
		if (fn_bb) t.push(bb[key][2], bb[key][3], bb[key][4], bb[key][5], bb[key][6], bb[key][7], bb[key][8], bb[key][9]);
		else t.push(len);
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
			if (c == '_inter') continue;
			//if (c == '_inter') continue;
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
		} else if (x['_inter'] != null && x['_inter'] >= len * 0.7) {
			type = 'inter';
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
}

function mg_cmd_paf2bl(args)
{
	var c, min_de = 0.01, max_de = 0.1, sub_de = 0.002, min_mapq = 5, min_len = 500, is_sub = false;
	while ((c = getopt(args, "d:s")) != null) {
		if (c == 'd') min_de = parseFloat(getopt.arg);
		else if (c == 's') is_sub = true;
	}
	if (args.length - getopt.ind < 1) {
		print("Usage: mgutils.js paf2bl <ins.paf>");
		print("Note: bedtk sub <(mgutils.js paf2bl ins.paf; cat bl100.bed) <(../mgutils.js paf2bl -s ins.paf) | bedtk merge");
		return;
	}
	var file = new File(args[getopt.ind]);
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var m, t = line.split("\t");
		if (/\ttp:A:[SI]/.test(line)) continue;
		if (parseInt(t[11]) < min_mapq) continue;
		if (parseInt(t[10]) < min_len) continue;
		if ((m = /\tde:f:(\S+)/.exec(line)) == null) continue;
		var de = parseFloat(m[1]);
		if (is_sub) {
			if (de > sub_de) continue;
		} else {
			if (de < min_de || de > max_de) continue;
		}
		print(t[5], t[7], t[8]);
		//print(line);
	}
	buf.destroy();
	file.close();
}

function mg_cmd_stableGaf(args)
{
	var c;
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length - getopt.ind < 1) {
		print("Usage: mgutils.js stableGaf <graph.gfa> <aln.gaf>");
		return;
	}

	var re = /\t(LN|SN|SO|SR):[Zi]:(\S+)/g;
	var file, buf = new Bytes();

	var pri_len = {}, segh = {};
	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^S\t(\S+)\t(\S+)(\t.*)/.exec(line)) == null) continue;
		var seg = m[1], len = m[2] == '*'? 0 : m[2].length, tags = m[3];
		var sn = null, so = -1, sr = -1;
		while ((m = re.exec(tags)) != null) {
			if (m[1] == "LN") len = parseInt(m[2]);
			else if (m[1] == "SN") sn = m[2];
			else if (m[1] == "SO") so = parseInt(m[2]);
			else if (m[1] == "SR") sr = parseInt(m[2]);
		}
		if (sn == null || so < 0 || sr < 0 || len <= 0)
			throw Error("failed to parse tags '" + tags + "'");
		segh[seg] = [sn, so, so + len, sr];
		if (sr == 0) {
			if (pri_len[sn] == null) pri_len[sn] = 0;
			pri_len[sn] = pri_len[sn] > so + len? pri_len[sn] : so + len;
		}
	}
	file.close();

	re = /([><])([^\s><]+)/g;
	file = args.length - getopt.ind < 2? new File() : new File(args[getopt.ind+1]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^(\S+)\t(\d+\t\d+\t\d+)\t([+-])\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t(.*)/.exec(line)) == null)
			continue;
		var s, a = [];
		while ((s = re.exec(m[4])) != null) {
			if (segh[s[2]] == null)
				throw Error("failed to find segment '" + s[2] + "'");
			var h = segh[s[2]], add_new = true;
			if (a.length) {
				var b = a[a.length - 1];
				if (b[0] == s[1] && h[3] == b[4] && h[0] == b[1]) {
					if (b[0] == '>') {
						if (h[1] == b[3]) b[3] = h[2], add_new = false;
					} else {
						if (h[2] == b[2]) b[2] = h[1], add_new = false;
					}
				}
			}
			if (add_new) a.push([s[1], h[0], h[1], h[2], h[3]]);
		}
		var path_len = 0, path = "";
		for (var i = 0; i < a.length; ++i)
			path_len += a[i][3] - a[i][2];
		if (path_len != parseInt(m[5]))
			throw Error("inconsistent path length for '" + m[1] + "': " + path_len + "!=" + m[5]);
		if (a.length == 1 && pri_len[a[0][1]] != null) {
			m[6] = parseInt(m[6]);
			m[7] = parseInt(m[7]);
			if (a[0][0] == '>') {
				m[6] += a[0][2], m[7] += a[0][2];
			} else {
				m[3] = m[3] == '+'? '-' : '+';
				var st = a[0][2] + (path_len - 1 - m[7]);
				var en = a[0][2] + (path_len - 1 - m[6]);
				m[6] = st, m[7] = en;
			}
			path_len = pri_len[a[0][1]];
			path = a[0][1];
		} else {
			var b = [];
			for (var i = 0; i < a.length; ++i)
				b.push(a[i][0] + a[i][1] + ':' + a[i][2] + '-' + a[i][3]);
			path = b.join("");
		}
		print(m[1], m[2], m[3], path, path_len, m[6], m[7], m[8]);
	}
	file.close();
	buf.destroy();
}

function mg_cmd_subgaf(args) // FIXME: this is BUGGY!!!
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
	var c, flank = 100, min_var_len = 100, min_test_len = 50, min_sc = 20.0, non_chr = false, out_err = false, flt_vcf = false;
	while ((c = getopt(args, "f:v:t:s:aeF")) != null) {
		if (c == 'f') flank = parseInt(getopt.arg);
		else if (c == 'v') min_var_len = parseInt(getopt.arg);
		else if (c == 't') min_test_len = parseInt(getopt.arg);
		else if (c == 's') min_sc = parseFloat(getopt.arg);
		else if (c == 'a') non_chr = true;
		else if (c == 'e') out_err = true;
		else if (c == 'F') flt_vcf = true;
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
		var flt = (t[6] != '.' && t[6] != 'PASS');
		if (flt_vcf && flt) continue;
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
			if (gt[i] == '.') continue;
			var x = parseInt(gt[i]);
			var l = al[x].length - ref.length;
			var x = l > 0? l : -l;
			if (max_diff < x) max_diff = x, max_ev = l;
		}
		if (max_diff < min_test_len) continue;
		if (vcf[t[0]] == null) vcf[t[0]] = [];
		vcf[t[0]].push([st, en, -1, max_diff, max_ev, flt, s[0]]);
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
		var ref_len = t[7] == '*'? 0 : t[7].length;
		var max_diff = 0, max_ev = 0;
		for (var i = 8; i < t.length; ++i) {
			var alt_len = t[i] == '*'? 0 : t[8].length;
			var l = alt_len - ref_len;
			var x = l > 0? l : -l;
			if (max_diff < x) max_diff = x, max_ev = l;
		}
		var st = parseInt(t[1]), en = parseInt(t[2]);
		rst[t[0]].push([st, en, -1, max_diff, max_ev]);
	}
	file.close();
	for (var ctg in rst) it_index(rst[ctg]);

	// sensitivity
	var n_vcf = [0, 0, 0], fn = [0, 0, 0];
	for (var ctg in vcf) {
		for (var i = 0; i < vcf[ctg].length; ++i) {
			var v = vcf[ctg][i];
			if (v[3] < min_var_len) continue;
			if (v[5]) continue;
			var st = v[0] - flank, en = v[1] + flank;
			if (st < 0) st = 0;
			if (!it_contained(bed[ctg], st, en)) continue;
			var sub = v[4] < 0? 1 : 2;
			++n_vcf[0], ++n_vcf[sub];
			var b = it_overlap(rst[ctg], st, en);
			if (b.length == 0) {
				if (out_err) print("FN", ctg, v[0], v[1], v[4], v[6]);
				++fn[0], ++fn[sub];
			}
		}
	}

	// specificity
	var n_rst = [0, 0, 0], fp = [0, 0, 0];
	for (var ctg in rst) {
		for (var i = 0; i < rst[ctg].length; ++i) {
			var v = rst[ctg][i];
			if (v[3] < min_var_len) continue;
			var st = v[0] - flank, en = v[1] + flank;
			if (st < 0) st = 0;
			if (!it_contained(bed[ctg], st, en)) continue;
			var sub = v[4] < 0? 1 : 2;
			++n_rst[0], ++n_rst[sub];
			var b = it_overlap(vcf[ctg], st, en);
			if (b.length == 0) {
				if (out_err) print("FP", ctg, v[0], v[1], v[4]);
				++fp[0], ++fp[sub];
			}
		}
	}

	print("NA", fn[0], n_vcf[0], (fn[0]/n_vcf[0]).toFixed(4));
	print("ND", fn[1], n_vcf[1], (fn[1]/n_vcf[1]).toFixed(4));
	print("NI", fn[2], n_vcf[2], (fn[2]/n_vcf[2]).toFixed(4));
	print("PA", fp[0], n_rst[0], (fp[0]/n_rst[0]).toFixed(4));
	print("PD", fp[1], n_rst[1], (fp[1]/n_rst[1]).toFixed(4));
	print("PI", fp[2], n_rst[2], (fp[2]/n_rst[2]).toFixed(4));
}

function mg_cmd_extractseg(args)
{
	function process(ctg, first, last, is_end) {
		if (ctg == null || first[0] == null || first[1] == null) return;
		if (first[0][7] == first[1][7]) return;
		if (first[0][7] < first[1][7]) {
			if (last[0][7] >= first[1][7]) return;
			if (is_end) print(ctg, last[0][8], first[1][7], '*', 0, '+');
			else print(ctg, last[0][7], first[1][8], '*', 0, '+');
		} else {
			if (last[1][7] >= first[0][7]) return;
			if (is_end) print(ctg, last[1][8], first[0][7], '*', 0, '-');
			else print(ctg, last[1][7], first[0][8], '*', 0, '-');
		}
	}

	var c, min_len = 100000, is_end = false;
	while ((c = getopt(args, "el:")) != null) {
		if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'e') is_end = true;
	}
	if (args.length - getopt.ind < 3) {
		print("Usage: mgutils.js extractseg <seg1> <seg2> <in.gaf> [...]");
		return;
	}

	var seg = [args[getopt.ind], args[getopt.ind+1]];
	var buf = new Bytes();
	for (var i = getopt.ind + 2; i < args.length; ++i) {
		var file = new File(args[i]);
		var flt = false;
		var first = [null, null], last = [null, null], ctg = null;
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (t[0] != "*") {
				process(ctg, first, last, is_end);
				flt = (parseInt(t[3]) - parseInt(t[2]) < min_len || parseInt(t[8]) - parseInt(t[7]) < min_len);
				first = [null, null];
				last = [null, null];
				ctg = t[0];
			} else if (!flt) {
				var s = t[1].substr(1);
				t[7] = parseInt(t[7]), t[8] = parseInt(t[8]);
				if (s == seg[0] && t[3] != '0') {
					if (first[0] == null) first[0] = t.slice(0);
					last[0] = t.slice(0);
				} else if (s == seg[1] && t[3] != '0') {
					if (first[1] == null) first[1] = t.slice(0);
					last[1] = t.slice(0);
				}
			}
		}
		process(ctg, first, last, is_end);
		file.close();
	}
	buf.destroy();
}

function mg_cmd_bed2sql(args)
{
	var c;
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length - getopt.ind == 0) {
		print("Usage: paste *.bed | mgutils.js bed2sql <sample.list> | sqlite3 rGFA.db");
		return;
	}

	var file, buf = new Bytes();

	var sample = [];
	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		sample.push(t[0]);
	}
	file.close();

	file = args.length - getopt.ind >= 2 && args[getopt.ind+1] != "-"? new File(args[getopt.ind+1]) : new File();
	print("DROP INDEX IF EXISTS idx_bwalk;");
	print("DROP INDEX IF EXISTS idx_cst;");
	print("DROP INDEX IF EXISTS idx_cen;");
	print("BEGIN TRANSACTION;");
	var wid = 0, bid = 0, ins_walk = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t.length != sample.length * 6)
			throw Error("Different number of samples");
		var h = {}, w = [], j = 0;
		for (var i = 5; i < t.length; i += 6, ++j) {
			if (t[i] == ".") continue;
			var s = t[i].split(":");
			if (!(s[0] in h)) {
				h[s[0]] = w.length;
				ins_walk.push([wid, bid, s[1], s[0]]);
				w.push([s[0], s[1], wid++]);
			}
			var v = [], x = w[h[s[0]]];
			v.push("'" + bid + "'", "'" + sample[j] + "'", "'" + x[2] + "'", "'" + s[3] + "'");
			v.push("'" + s[4] + "'", "'" + s[5] + "'", "'" + (s[2] == '+'? 1 : -1) + "'");
			print("INSERT INTO call (bid,sample,wid,ctg,start,end,strand) VALUES (" + v.join(",") + ");");
		}
		++bid;
	}
	for (var i = 0; i < ins_walk.length; ++i) {
		var w = ins_walk[i], v = [];
		for (var j = 0; j < w.length; ++j)
			v.push("'" + w[j] + "'");
		print("INSERT INTO bwalk (wid,bid,len,walk) VALUES (" + v.join(",") + ");");
	}
	print("END TRANSACTION;");
	print("CREATE INDEX IF NOT EXISTS idx_bwalk ON bwalk (bid);");
	print("CREATE INDEX IF NOT EXISTS idx_cst   ON call  (ctg, start);");
	print("CREATE INDEX IF NOT EXISTS idx_cen   ON call  (ctg, end);");
	file.close();

	buf.destroy();
}

function mg_cmd_merge(args)
{
	var c, fn_anno = null;
	while ((c = getopt(args, "a:")) != null) {
		if (c == 'a') fn_anno = getopt.arg;
	}
	if (args.length - getopt.ind == 0) {
		print("Usage: paste *.bed | mgutils.js merge [-a <anno.bed>] -");
		return;
	}

	var file, buf = new Bytes();
	var anno = {};
	if (fn_anno) {
		file = new File(fn_anno);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var key = [t[0], t[1], t[2]].join("_");
			anno[key] = t[11];
		}
		file.close();
	}
	file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	print('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">');
	print('##INFO=<ID=NA,Number=1,Type=Integer,Description="Number of alleles">');
	print('##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count">');
	print('##INFO=<ID=ALEN,Number=.,Type=Integer,Description="Length of each allele">');
	print('##INFO=<ID=ANNO,Number=1,Type=String,Description="Annotation">');
	print('##INFO=<ID=VS,Number=1,Type=String,Description="Start vertex">');
	print('##INFO=<ID=VE,Number=1,Type=String,Description="End vertex">');
	print('##INFO=<ID=AWALK,Number=.,Type=String,Description="Walk of each allele">');
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">');
	print('##FORMAT=<ID=CSTRAND,Number=1,Type=String,Description="Contig strand">');
	print('##FORMAT=<ID=CTG,Number=1,Type=String,Description="Contig name">');
	print('##FORMAT=<ID=CS,Number=1,Type=String,Description="Contig start, BED-like">');
	print('##FORMAT=<ID=CE,Number=1,Type=String,Description="Contig end, BED-like">');
	print("#CHROM\tSTART\tEND\tINFO\tFORMAT");
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var a = [t[0], t[1], t[2], "", "GT:CSTRAND:CTG:CS:CE"];
		var ah = {}, aa = [], b = [], ns = 0;
		for (var j = 5; j < t.length; j += 6) {
			if (t[j] == ".") {
				b.push(["."]);
				continue;
			}
			++ns;
			var s = t[j].split(":");
			if (ah[s[0]] == null) {
				ah[s[0]] = aa.length;
				aa.push({walk:s[0], len:s[1], cnt:0});
			}
			var k = ah[s[0]];
			++aa[k].cnt;
			s[0] = k;
			b.push(s);
		}
		for (var i = 0; i < aa.length; ++i)
			aa[i].i = i;
		aa.sort(function(a,b) { return b.cnt - a.cnt });
		var i2a = [], alen = [], awalk = [], ac = [];
		for (var i = 0; i < aa.length; ++i) {
			i2a[aa[i].i] = i;
			alen[i] = aa[i].len;
			awalk[i] = aa[i].walk;
			ac[i] = aa[i].cnt;
		}
		for (var j = 0; j < b.length; ++j) {
			if (b[j][0] != ".") {
				var i = b[j].shift();
				b[j][0] = i2a[i];
				a.push(b[j].join(":"));
			} else a.push(".");
		}
		var info = ["NS="+ns, "NA="+aa.length, "ALEN="+alen.join(","), "AC="+ac.join(",")];
		var key = [t[0], t[1], t[2]].join("_");
		if (anno[key] != null) info.push("ANNO="+anno[key]);
		info.push("VS="+t[3], "VE="+t[4], "AWALK="+awalk.join(","));
		a[3] = info.join(";");
		print(a.join("\t"));
	}
	buf.destroy();
	file.close();
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: mgutils.js <command> [arguments]");
		print("Commands:");
		print("  stableGaf    convert unstable GAF to stable GAF");
		print("  renamefa     add a prefix to sequence names in FASTA");
		print("  paf2bl       blacklist regions from insert-to-ref alignment");
		print("  anno         annotate short sequences");
		print("  extractseg   extract a segment from GAF");
		print("  merge        merge per-sample --call BED");
		print("  bed2sql      generate SQL from --call BED");
		//print("  subgaf       extract GAF overlapping with a region (BUGGY)");
		//print("  sveval       evaluate SV accuracy");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'renamefa') mg_cmd_renamefa(args);
	else if (cmd == 'paf2bl') mg_cmd_paf2bl(args);
	else if (cmd == 'anno') mg_cmd_anno(args);
	else if (cmd == 'subgaf') mg_cmd_subgaf(args);
	else if (cmd == 'sveval') mg_cmd_sveval(args);
	else if (cmd == 'joinfa') mg_cmd_joinfa(args);
	else if (cmd == 'stableGaf') mg_cmd_stableGaf(args);
	else if (cmd == 'bed2sql') mg_cmd_bed2sql(args);
	else if (cmd == 'extractseg') mg_cmd_extractseg(args);
	else if (cmd == 'merge') mg_cmd_merge(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
