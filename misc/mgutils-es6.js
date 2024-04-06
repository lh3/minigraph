#!/usr/bin/env k8

const version = "r578";

/**************
 * From k8.js *
 **************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

/***************
 * Subcommands *
 ***************/

function mg_cmd_merge2vcf(args) {
	let opt = { max_allele:15, ref_index:0, fn_sample:null, sample:[] };
	for (const o of getopt(args, "r:a:s:", [])) {
		if (o.opt == "-r") opt.ref_index = parseInt(o.arg);
		else if (o.opt == "-a") opt.max_allele = parseInt(o.arg);
		else if (o.opt == "-s") opt.fn_sample = o.arg;
	}
	if (args.length == 0) {
		print(`Usage: mgutils-es6.js merge2vcf [options] <in.bed>`);
		print(`Options:`);
		print(`  -r INT    which sample corresponds to the reference [${opt.ref_index}]`);
		print(`  -a INT    max allele number [${opt.max_allele}]`);
		print(`  -s FILE   list of sample names, one per line []`);
		return;
	}
	let file, buf = new Bytes();
	if (opt.fn_sample) {
		file = new File(opt.fn_sample);
		while (file.readline(buf) >= 0) {
			const t = buf.toString().split(/\s+/);
			opt.sample.push(t[0]);
		}
		file.close();
	}
	file = new File(args[0]);
	let hdr = [];
	hdr.push(`##fileformat=VCFv4.2`);
	hdr.push(`##ALT=<ID=CNV,Description="description">`);
	hdr.push(`##FORMAT=<ID=GT0,Number=1,Type=String,Description="Original genotype">`);
	for (let i = 1; i <= opt.max_allele; ++i)
		hdr.push(`##ALT=<ID=X:${i},Description="Allele ${i}">`);
	let n_sample = opt.sample.length;
	while (file.readline(buf) >= 0) {
		let line = buf.toString();
		if (line[0] == "#" && line[1] == "#") {
			hdr.push(line);
		} else if (line[0] == '#') {
			let t = line.split("\t");
			let a = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"];
			if (t.length <= 5) {
				if (opt.sample.length == 0)
					throw Error("No samples found. Please provide sample names with option '-s'");
				for (let i = 0; i < opt.sample.length; ++i)
					a.push(opt.sample[i]);
			} else {
				for (let i = 6; i < t.length; ++i)
					a.push(t[i]);
			}
			for (let i = 0; i < hdr.length; ++i)
				print(hdr[i]);
			print(`#${a.join("\t")}`);
		} else {
			let t = buf.toString().split("\t");
			if (n_sample == 0) n_sample = t.length - 5;
			if (n_sample != t.length - 5) throw Error("different number of samples");
			let a = [t[0], t[1], ".", "N", "", "30", "PASS"]
			let m, ref = -1;
			if ((m = /^(\d+)/.exec(t[5 + opt.ref_index])) != null)
				ref = parseInt(m[1]);
			if ((m = /\bNA=(\d+)/.exec(t[3])) == null) throw Error("No NA tag");
			let na = parseInt(m[1]), a2v = [];
			for (let i = 0; i < na; ++i)
				a2v[i] = i;
			if (ref >= 0) {
				for (let i = 0; i < ref; ++i)
					a2v[i] = i + 1;
				a2v[ref] = 0;
			}
			let al = [];
			for (let i = 1; i < na && i <= opt.max_allele; ++i)
				al.push(`<X:${i}>`);
			a[4] = al.length? al.join(",") : ".";
			let info = [`END=${t[2]}`];
			const re = /([^\s=;]+)=([^\s=;]+)/g;
			while ((m = re.exec(t[3])) != null) {
				if (m[1] == "ALEN" || m[1] == "AWALK" || m[1] == "AC") {
					const s = m[2].split(",");
					if (s.length != na) throw Error("Inconsistent number of alleles");
					let p = [];
					if (m[1] == "AC") {
						for (let i = 0; i < s.length; ++i)
							if (a2v[i] != 0)
								p.push(s[i]);
					} else {
						for (let i = 0; i < s.length; ++i)
							p[a2v[i]] = s[i];
					}
					if (m[1] != "AC" || p.length > 0)
						info.push(`${m[1]}=${p.join(",")}`);
				} else if (m[1] == "NS") {
					info.push(`AN=${m[2]}`);
					info.push(`${m[1]}=${m[2]}`);
				} else {
					info.push(`${m[1]}=${m[2]}`);
				}
			}
			a.push(info.join(";"), "GT:GT0");
			for (let i = 5; i < t.length; ++i) {
				if (t[i] == ".") {
					a.push(".");
				} else if ((m = /^(\d+)(\S*)/.exec(t[i])) != null) {
					const al = a2v[parseInt(m[1])];
					const al_cap = al < opt.max_allele? al : opt.max_allele;
					a.push(`${al_cap}:${al}`);
				}
			}
			print(a.join("\t"));
		}
	}
	file.close();
	buf.destroy();
}

function mg_cmd_addsample(args) {
	if (args.length < 2) {
		print("Usage: mgutils-es6.js addsample <merged.bed> <sample.txt>");
		return;
	}
	let file, buf = new Bytes();
	file = new File(args[1]);
	let sample = [];
	while (file.readline(buf) >= 0) {
		let t = buf.toString().split(/\s+/);
		sample.push(t[0]);
	}
	file.close();
	file = new File(args[0]);
	while (file.readline(buf) >= 0) {
		const line = buf.toString();
		if (line[0] != "#" || (line[0] == "#" && line[1] == "#")) {
			print(buf);
		} else {
			print("#CHROM", "START", "END", "INFO", "FORMAT", sample.join("\t"));
		}
	}
	file.close();
	buf.destroy();
}

function mg_revcomp(s) {
	function complement(x) { return { a:'t', t:'a', g:'c', c:'g' }[x] }
	return s.split('').reverse().map(complement).join('');
}

function mg_cmd_getsv(args) {
	let opt = { min_mapq:5, min_mapq_end:30, min_frac:0.7, min_len:100, min_aln_len_end:2000, min_aln_len_mid:50, max_cnt:5, dbg:false, polyA_pen:5, polyA_drop:100, name:"foo", cen:{} };
	for (const o of getopt(args, "q:Q:l:dc:a:e:m:n:b:", [])) {
		if (o.opt == "-q") opt.min_mapq = parseInt(o.arg);
		else if (o.opt == "-Q") opt.min_mapq_end = parseInt(o.arg);
		else if (o.opt == "-l") opt.min_len = parseInt(o.arg);
		else if (o.opt == "-d") opt.dbg = true;
		else if (o.opt == "-f") opt.min_frac = parseFloat(o.arg);
		else if (o.opt == "-c") opt.max_cnt = parseInt(o.arg);
		else if (o.opt == "-a") opt.polyA_pen = parseInt(o.arg);
		else if (o.opt == "-e") opt.min_aln_len_end = parseInt(o.arg);
		else if (o.opt == "-m") opt.min_aln_len_mid = parseInt(o.arg);
		else if (o.opt == "-n") opt.name = o.arg;
		else if (o.opt == "-b") {
			for (const line of k8_readline(o.arg)) {
				const t = line.split("\t");
				if (opt.cen[t[0]] == null) opt.cen[t[0]] = [];
				opt.cen[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
			}
			for (const ctg in opt.cen)
				opt.cen[ctg].sort(function(a,b) { return a[0]-b[0] });
		}
	}
	if (args.length == 0) {
		print("Usage: mgutils-es6.js getsv [options] <stable.gaf>");
		print("Options:");
		print(`  -n STR     sample name [${opt.name}]`);
		print(`  -q INT     min mapq [${opt.min_mapq}]`);
		print(`  -l INT     min INDEL len [${opt.min_len}]`);
		print(`  -f FLOAT   min mapped query fraction [${opt.min_frac}]`);
		print(`  -c INT     max number of long INDELs per read [${opt.max_cnt}]`);
		print(`  -a INT     penalty for non-polyA bases [${opt.polyA_pen}]`);
		print(`  -Q INT     min mapq for alignment ends [${opt.min_mapq_end}]`);
		print(`  -e INT     min alignment length at ends [${opt.min_aln_len_end}]`);
		print(`  -m INT     min alignment length in the middle [${opt.min_aln_len_mid}]`);
		print(`  -b FILE    BED for centromeres []`);
		return;
	}

	let re = /(\d+)([=XIDMSHN])/g; // regex for CIGAR
	let re_path = /([><])([^><:\s]+):(\d+)-(\d+)/g; // regex for path/ctg
	let re_ds = /([\+\-\*:])([A-Za-z\[\]0-9]+)/g; // regex for the ds tag
	let re_tsd = /(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?/; // regex for parsing TSD

	function cal_cen_dist(opt, ctg, pos) {
		if (opt.cen[ctg] == null) return 1e9;
		let min = 1e9;
		for (let i = 0; i < opt.cen[ctg].length; ++i) { // TODO: binary search would be better
			const b = opt.cen[ctg][i];
			const d = pos < b[0]? b[0] - pos : pos < b[1]? 0 : pos - b[1];
			min = min < d? min : d;
		}
		return min;
	}

	function cal_cen_overlap(opt, ctg, st0, en0) {
		if (opt.cen[ctg] == null) return 0;
		let cov_st = 0, cov_en = 0, cov = 0;
		for (let i = 0; i < opt.cen[ctg].length; ++i) { // TODO: binary search would be better
			const b = opt.cen[ctg][i];
			if (b[1] <= st0 || b[0] >= en0) continue; // not overlapping with [st0, en0)
			const st1 = b[0] > st0? b[0] : st0;
			const en1 = b[1] < en0? b[1] : en0;
			if (st1 > cov_en) {
				cov += cov_en - cov_st;
				cov_st = st1, cov_en = en1;
			} else cov_en = cov_en > en1? cov_en : en1;
		}
		cov += cov_en - cov_st;
		return cov;
	}

	/*******************************************
	 * Extract long indels contained in CIGARs *
	 *******************************************/

	function cal_polyA_len(opt, int_seq) {
		let polyA_len = 0, polyT_len = 0, polyA_max = 0, polyT_max = 0;
		let score, max, max_j;
		// look for polyA on the 3'-end
		score = max = 0, max_j = int_seq.length;
		for (let j = int_seq.length - 1; j >= 0; --j) {
			if (int_seq[j] == 'A' || int_seq[j] == 'a') ++score;
			else score -= opt.polyA_pen;
			if (score > max) max = score, max_j = j;
			else if (max - score > opt.polyA_drop) break;
		}
		polyA_len = int_seq.length - max_j, polyA_max = max;
		// look for polyT on the 5'-end
		score = max = 0, max_j = -1;
		for (let j = 0; j < int_seq.length; ++j) {
			if (int_seq[j] == 'T' || int_seq[j] == 't') ++score;
			else score -= opt.polyA_pen;
			if (score > max) max = score, max_j = j;
			else if (max - score > opt.polyA_drop) break;
		}
		polyT_len = max_j + 1, polyT_max = max;
		// choose the longer one
		return polyA_max >= polyT_max? polyA_len : -polyT_len;
	}

	function get_indel(opt, z) {
		if (z.length == 0) return;
		for (let j = 0; j < z.length; ++j) {
			const y = z[j];
			if (y.qen - y.qst < y.qlen * opt.min_frac) continue; // ignore short alignments
			let m, a = [], x = y.tst;
			while ((m = re.exec(y.cg)) != null) { // collect the list of long indels
				const op = m[2], len = parseInt(m[1]);
				if (len >= opt.min_len) {
					if (op === "I")
						a.push({ st:x-1, en:x+1,   len:len,  indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:"." });
					else if (op === "D")
						a.push({ st:x,   en:x+len, len:-len, indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:"." });
				}
				if (op == "M" || op == "=" || op == "X" || op == "D")
					x += len;
			}
			if (a.length == 0 || a.length > opt.max_cnt) continue;
			// parse ds:Z
			if (y.ds) { // this MUST match CIGAR parsing
				let i = 0, x = y.tst, m;
				while ((m = re_ds.exec(y.ds)) != null) {
					const op = m[1], str = m[2];
					const seq = op === "+" || op === "-"? str.replace(/[\[\]]/g, "") : "";
					const len = op === ":"? parseInt(str) : op === "*"? 1 : op === "+" || op === "-"? seq.length : -1;
					if (len < 0) throw Error("can't determine length from the ds tag");
					if (len >= opt.min_len) { // extract INDEL sequence and check consistency with CIGAR
						if (op === "+") {
							if (a[i].st != x - 1 || a[i].en != x + 1 || a[i].len != len)
								throw Error(`inconsistent insertion at line ${lineno}`);
							a[i++].indel_seq = str;
						} else if (op === "-") {
							if (a[i].st != x || a[i].en != x + len || a[i].len != -len)
								throw Error(`inconsistent deletion at line ${lineno}`);
							a[i++].indel_seq = str;
						}
					}
					if (op === "*" || op === ":" || op === "-")
						x += len;
				}
				for (let i = 0; i < a.length; ++i) { // compute TSD and polyA lengths
					if ((m = re_tsd.exec(a[i].indel_seq)) == null)
						throw Error("Bug!");
					const tsd = (m[5]? m[5] : "") + (m[2]? m[2] : "");
					const int_seq = m[3]; // internal sequence
					a[i].tsd_len = tsd.length;
					a[i].tsd_seq = tsd;
					a[i].int_seq = int_seq;
					if (int_seq.length > 0)
						a[i].polyA_len = cal_polyA_len(opt, int_seq);
				}
			} // ~if(y.ds)
			if (opt.dbg) print('X0', line);
			let seg = []; // reference segments in the path
			if (/[><]/.test(y.path)) { // with ><: this is a path
				let x = 0;
				if (y.strand != '+') throw Error("reverse strand on path");
				while ((m = re_path.exec(y.path)) != null) {
					const st = parseInt(m[3]), en = parseInt(m[4]);
					seg.push([m[2], st, en, m[1] == '>'? 1 : -1, x, x + (en - st)]);
					x += en - st;
				}
			} else { // this is a contig name
				seg.push([y.path, 0, y.tlen, 1, 0, y.tlen]);
			}
			let st = [], en = [];
			for (let i = 0, k = 0; i < a.length; ++i) { // start
				while (k < seg.length && seg[k][5] <= a[i].st) ++k;
				if (k == seg.length) throw Error("failed to find start");
				const l = a[i].st - seg[k][4];
				if (seg[k][3] > 0)
					st.push([k, seg[k][1] + l, seg[k][1] + l]);
				else
					st.push([k, seg[k][2] - l, seg[k][1] + l]);
			}
			for (let i = 0, k = 0; i < a.length; ++i) { // end
				while (k < seg.length && seg[k][5] <= a[i].en) ++k;
				if (k == seg.length) throw Error("failed to find end");
				const l = a[i].en - seg[k][4];
				if (seg[k][3] > 0)
					en.push([k, seg[k][1] + l, seg[k][1] + l]);
				else
					en.push([k, seg[k][2] - l, seg[k][1] + l]);
			}
			for (let i = 0; i < a.length; ++i) {
				if (opt.dbg) print('X2', a[i].st, a[i].en, st[i][0], en[i][0]);
				if (st[i][0] === en[i][0] && seg[st[i][0]][3] < 0) { // then reverse complement tsd, polyA and insert
					a[i].polyA_len = -a[i].polyA_len;
					a[i].tsd_seq = mg_revcomp(a[i].tsd_seq);
					a[i].int_seq = mg_revcomp(a[i].int_seq);
				}
				let info1 = (a[i].len > 0? "SVTYPE=INS" : "SVTYPE=DEL") + `;SVLEN=${a[i].len};tsd_len=${a[i].tsd_len};polyA_len=${a[i].polyA_len}`;
				const info2 = `source=${opt.name};tsd_seq=${a[i].tsd_seq.length>0?a[i].tsd_seq:"."};insert=${a[i].int_seq.length>0?a[i].int_seq:"."}`;
				if (st[i][0] == en[i][0]) { // on the same segment
					const s = seg[st[i][0]];
					const strand2 = s[3] > 0? y.strand : y.strand === '+'? '-' : '+';
					if (opt.cen[s[0]] != null) {
						const dist_st = cal_cen_dist(opt, s[0], st[i][1]);
						const dist_en = cal_cen_dist(opt, s[0], en[i][1]);
						info1 += `;cen_dist=${dist_st < dist_en? dist_st : dist_en}`
					}
					print(s[0], st[i][1] < en[i][1]? st[i][1] : en[i][1], st[i][1] > en[i][1]? st[i][1] : en[i][1], y.qname, y.mapq, strand2, `${info1};${info2}`);
				} else { // on different segments
					let path = [], len = 0;
					for (let j = st[i][0]; j <= en[i][0]; ++j) {
						const s = seg[j];
						len += s[2] - s[1];
						path.push((s[3] > 0? '>' : '<') + s[0] + `:${s[1]}-${s[2]}`);
					}
					const off = seg[st[i][0]][4];
					print(path.join(""), a[i].st - off, a[i].en - off, y.qname, y.mapq, '+', `${info1};${info2}`);
				}
			} // ~for(i)
		} // ~for(j)
	} // ~get_indel()

	/*********************************
	 * Extract alignment breakpoints *
	 *********************************/

	function get_end_coor(y) {
		let r = [{}, {}];
		if (/^[><]/.test(y.path)) { // a path
			if (y.strand != '+') throw Error("reverse strand on path");
			let x = 0, i = 0;
			while ((m = re_path.exec(y.path)) != null) {
				const st = parseInt(m[3]), en = parseInt(m[4]), len = en - st;
				if (y.tst >= x && y.tst < x + len) {
					r[0] = { ctg:m[2], ori:m[1], pos:-1 };
					r[0].pos = m[1] === ">"? st + (y.tst - x) : st + (x + len - y.tst) - 1;
				}
				if (y.ten > x && y.ten <= x + len) {
					r[1] = { ctg:m[2], ori:m[1], pos:-1 };
					r[1].pos = m[1] === ">"? st + (y.ten - x) - 1 : st + (x + len - y.ten);
				}
				x += len;
			}
		} else { // a contig
			r[0] = { ctg:y.path, ori: y.strand === "+"? ">" : "<", pos:-1 };
			r[0].pos = y.strand === "+"? y.tst : y.ten - 1;
			r[1] = { ctg:y.path, ori: y.strand === "+"? ">" : "<", pos:-1 };
			r[1].pos = y.strand === "+"? y.ten - 1 : y.tst;
		}
		r[0].ql = r[1].ql = y.qen - y.qst;
		return r;
	}

	function infer_svtype(opt, c0, c1, ori, qgap) { // NB: c0 MUST have the smaller coordinate
		if (c0.ctg != c1.ctg) return { st:-1, en:-1, str:"SVTYPE=BND" };
		const l = c1.pos - c0.pos + 1;
		if (l < 0) throw Error("Bug!");
		if (ori === ">>" && qgap < l && l - qgap >= opt.min_len) { // deletion
			const st = qgap < 0? c0.pos + qgap : c0.pos;
			const en = qgap < 0? c1.pos + 1 - qgap : c1.pos + 1;
			return { st:st, en:en, str:`SVTYPE=DEL;SVLEN=${-(l - qgap)};sv_region=${st},${en};tsd_len=${qgap < 0? -qgap : 0}` };
		}
		if (ori === ">>" && l < qgap && qgap - l >= opt.min_len) // insertion without TSD
			return { st:c0.pos, en:c1.pos+1, str:`SVTYPE=INS;SVLEN=${qgap - l};sv_region=${c0.pos},${c1.pos+1}` };
		if (ori === "<<" && qgap > 0 && (l < c0.ql || l < c1.ql) && qgap + l >= opt.min_len) // insertion with TSD
			return { st:c0.pos, en:c1.pos+1, str:`SVTYPE=INS;SVLEN=${qgap + l};sv_region=${c0.pos},${c1.pos+1};tsd_len=${l}` }; // TODO: is sv_region correct?
		if (ori === "<<" && qgap + l >= opt.min_len) { // tandem duplication; similar to insertion with TSD
			const st = qgap < 0? c0.pos : c0.pos > qgap? c0.pos - qgap : 0;
			const en = qgap < 0? c1.pos + 1 : c1.pos + 1 + qgap;
			return { st:st, en:en, str:`SVTYPE=DUP;SVLEN=${qgap + l};sv_region=${st},${en}` };
		}
		if ((ori === "<>" || ori === "><") && l >= opt.min_len) { // inversion
			const st = qgap < 0? c0.pos + qgap : c0.pos;
			const en = qgap < 0? c1.pos + 1 - qgap : c1.pos + 1;
			return { st:st, en:en, str:`SVTYPE=INV;SVLEN=${l - qgap};sv_region=${st},${en}` };
		}
		return { st:-1, en:-1, str:"SVTYPE=BND" };
	}

	function get_breakpoint(opt, z) {
		if (z.length < 2) return;
		z.sort(function(a,b) { return a.qst - b.qst }); // sort by start position on the read
		// filter out short alignment towards the end of the read
		let zen = z.length;
		for (let j = z.length - 1; j >= 0; --j) {
			const y = z[j];
			if (y.qen - y.qst < opt.min_aln_len_end || y.mapq < opt.min_mapq_end) zen = j;
			else break;
		}
		if (zen < 2) return;
		// filter out short alignment towards the start of the read
		let zst = 0;
		for (let j = 0; j < zen; ++j) {
			const y = z[j];
			if (y.qen - y.qst < opt.min_aln_len_end || y.mapq < opt.min_mapq_end) zst = j + 1;
			else break;
		}
		if (zen - zst < 2) return;
		// construct the final alignment list
		let zz = [];
		for (let j = zst; j < zen; ++j)
			if (z[j].qen - z[j].qst >= opt.min_aln_len_mid)
				zz.push(z[j]);
		if (zz.length < 2) return; // shouldn't happen if mid<end; just in case
		// calculate the end coordinates
		for (let j = 0; j < zz.length; ++j) {
			const r = get_end_coor(zz[j]);
			zz[j].coor = r;
		}
		// extract alignment breakpoints
		for (let j = 1; j < zz.length; ++j) {
			const y0 = zz[j-1], y1 = zz[j];
			const qgap = y1.qst - y0.qen;
			let c0 = y0.coor[1], c1 = y1.coor[0], strand2 = "+", ori = c0.ori + c1.ori;
			if (!(c0.ctg < c1.ctg || (c0.ctg === c1.ctg && c0.pos < c1.pos)))
				c0 = y1.coor[0], c1 = y0.coor[1], strand2 = "-", ori = (c1.ori === ">"? "<" : ">") + (c0.ori === ">"? "<" : ">");
			const sv_info = infer_svtype(opt, c0, c1, ori, qgap);
			let cen_str = "";
			if (opt.cen[c0.ctg] != null || opt.cen[c1.ctg] != null) {
				const dist0 = cal_cen_dist(opt, c0.ctg, c0.pos);
				const dist1 = cal_cen_dist(opt, c1.ctg, c1.pos);
				cen_str = `;cen_dist=${dist0<dist1?dist0:dist1}`;
				if (sv_info.st >= 0 && sv_info.en >= sv_info.st) {
					const ov = cal_cen_overlap(opt, c0.ctg, sv_info.st, sv_info.en);
					cen_str += `;cen_overlap=${ov}`;
				}
			}
			print(c0.ctg, c0.pos, ori, c1.ctg, c1.pos, y0.qname, y0.mapq < y1.mapq? y0.mapq : y1.mapq, strand2,
				  `${sv_info.str};qgap=${qgap};mapq=${y0.mapq},${y1.mapq};aln_len=${y0.qen-y0.qst},${y1.qen-y1.qst}${cen_str};source=${opt.name}`);
		}
	} // ~get_breakpoint()

	let lineno = 0, z = [];
	for (const line of k8_readline(args[0])) {
		++lineno;
		let t = line.split("\t");
		if (t.length < 11) continue; // SAM has 11 columns at least; PAF has 12 columns at least
		if (z.length > 0 && t[0] != z[0].qname) {
			get_indel(opt, z);
			get_breakpoint(opt, z);
			z = [];
		}
		// parse format
		let y = { qname:t[0], mapq:0, qst:-1, qen:-1, qlen:-1, tlen:-1, tst:-1, cg:null, ds:null, path:null, strand:null };
		if (t.length >= 12 && (t[4] === "+" || t[4] === "-")) { // parse PAF or GAF
			y.mapq = parseInt(t[11]);
			if (y.mapq < opt.min_mapq) continue;
			y.qlen = parseInt(t[1]);
			y.qst = parseInt(t[2]);
			y.qen = parseInt(t[3]);
			y.strand = t[4];
			y.path = t[5];
			y.tlen = parseInt(t[6]);
			y.tst = parseInt(t[7]);
			y.ten = parseInt(t[8]);
			let tp = "";
			for (let i = 12; i < t.length; ++i) {
				if (t[i].substr(0, 5) === "cg:Z:")
					y.cg = t[i].substr(5);
				else if (t[i].substr(0, 5) === "ds:Z:")
					y.ds = t[i].substr(5);
				else if (t[i].substr(0, 5) === "tp:A:")
					tp = t[i].substr(5);
			}
			if (tp != "P") continue; // filter out secondary alignments
			if (y.cg == null) continue;
		} else { // parse SAM
			if (t[0][0] === "@") continue;
			const flag = parseInt(t[1]);
			if (flag & 0x100) continue;
			y.mapq = parseInt(t[4]);
			if (y.mapq < opt.min_mapq) continue;
			y.strand = (flag&0x10)? "-" : "+";
			y.path = t[2];
			y.tlen = 0xffffffff; // tlen doesn't need to be accurate for SAM or PAF
			y.tst = parseInt(t[3]) - 1;
			y.ten = -1;
			y.cg = t[5];
			let m;
			y.qst = (m = /^(\d+)[SH]/.exec(cg)) != null? parseInt(m[1]) : 0;
			y.qlen = 0;
			while ((m = re.exec(cg)) != null) {
				const op = m[2];
				if (op == "S" || op == "H" || op == "M" || op == "=" || op == "I")
					y.qlen += parseInt(m[1]);
			}
			for (let i = 11; i < t.length; ++i)
				if (t[i].substr(0, 5) === "ds:Z:")
					y.ds = t[i].substr(5);
			y.qen = y.qlen - ((m = /(\d+)[SH]$/.exec(cg)) != null? parseInt(m[1]) : 0);
		}
		z.push(y);
	}
	get_indel(opt, z);
	get_breakpoint(opt, z);
}

function mg_cmd_mergesv(args) {
	let opt = { min_cnt:3, min_cnt_rt:1, min_rt_len:10, win_size:100, max_diff:0.05, min_cen_dist:500000 };
	for (const o of getopt(args, "w:d:c:e:r:")) {
		if (o.opt === "-w") opt.win_size = parseInt(o.arg);
		else if (o.opt === "-d") opt.max_diff = parseInt(o.arg);
		else if (o.opt === "-c") opt.min_cnt = parseInt(o.arg);
		else if (o.opt === "-e") opt.min_cen_dist = parseInt(o.arg);
		else if (o.opt === "-r") opt.min_rt_len = parseInt(o.arg);
	}
	if (args.length == 0) {
		print("Usage: sort -k1,1 -k2,2n getsv.txt | mgutils-es6.js mergesv [options] -");
		print("Options:");
		print(`  -c INT     min read count [${opt.min_cnt}]`);
		print(`  -w INT     window size [${opt.win_size}]`);
		print(`  -d FLOAT   max allele length difference ratio [${opt.max_diff}]`);
		print(`  -e INT     min distance to centromeres [${opt.min_cen_dist}]`);
		print(`  -r INT     min min(TSD_len,polyA_len) to tag a candidate RT [${opt.min_rt_len}]`);
		return;
	}

	function parse_sv(t) {
		const re_info = /([^;\s=]+)=([^;\s=]+)/g;
		let v = { ctg:t[0], pos:parseInt(t[1]), st:-1, en:-1, ori:null, ctg2:null, pos2:null, _mapq:0, strand:null, is_bp:false, info:null }
		v.is_bp = /[><]/.test(t[2]);
		const off = v.is_bp? 6 : 4; // offset of the mapq column
		v._mapq = parseInt(t[off]);
		v.strand = t[off+1];
		v.info = t[off+2];
		let m;
		while ((m = re_info.exec(t[off+2])) != null) // parse INFO
			v[m[1]] = m[2];
		if (v.SVTYPE == null || v.source == null)
			throw Error("missing SVTYPE or source");
		if (v.SVLEN != null) v.SVLEN = parseInt(v.SVLEN);
		if (v.cen_dist != null) v.cen_dist = parseInt(v.cen_dist);
		if (v.cen_overlap != null) v.cen_overlap = parseInt(v.cen_overlap);
		if (!v.is_bp) {
			v.st = v.pos, v.en = parseInt(t[2]);
		} else {
			v.ori = t[2], v.ctg2 = t[3], v.pos2 = parseInt(t[4]);
			if (v.sv_region != null && (m = /(\d+),(\d+)/.exec(v.sv_region)) != null)
				v.st = parseInt(m[1]), v.en = parseInt(m[2]);
		}
		return v;
	}

	function same_sv(opt, v, w) {
		if (v.is_bp != w.is_bp) return false;
		if (v.ctg != w.ctg) return false; // not on the same contig
		if (v.SVTYPE != w.SVTYPE) return false; // not the same type
		if (v.is_bp && w.is_bp && v.ori != w.ori) { // test inversions
			if (!((v.ori === "><" && w.ori === "<>") || (v.ori === "<>" || w.ori === "><")))
				return false;
		}
		if (v.pos - w.pos > opt.win_size || w.pos - v.pos > opt.win_size) return false; // pos differ too much
		if (!v.is_bp) {
			if (v.en - w.en > opt.win_size || w.en - v.en > opt.win_size) return false; // end differ too much
		} else {
			if (v.ctg2 != w.ctg2) return false;
			if (v.pos2 - w.pos2 > opt.win_size || w.pos2 - v.pos2 > opt.win_size) return false;
		}
		if (v.SVLEN != null && w.SVLEN != null) {
			if (v.SVLEN * w.SVLEN <= 0) return false; // redundant but doesn't hurt to check
			const vl = Math.abs(v.SVLEN);
			const wl = Math.abs(w.SVLEN);
			if (Math.abs(vl - wl) > .5 * (vl + wl) * opt.max_diff) return false; // SVLEN differ too much
		}
		return true; // TODO: probably we don't want to check tsd_len as it may be cut short by a sequencing error
	}

	function write_sv(opt, s) {
		if (s.length == 0) return;
		const v = s[s.length>>1];
		// filter by cen_dist
		if (opt.min_cen_dist > 0) {
			if (v.cen_overlap != null && v.cen_overlap > 0) return;
			if (v.cen_dist != null && v.cen_dist <= opt.min_cen_dist) return;
		}
		// calculate rt_len
		let rt_len_arr = [], rt_len = 0;
		for (let i = 0; i < s.length; ++i)
			if (s[i].tsd_len != null && s[i].polyA_len != null)
				rt_len_arr.push(s[i].tsd_len < Math.abs(s[i].polyA_len)? s[i].tsd_len : Math.abs(s[i].polyA_len));
		if (rt_len_arr.length > 0)
			rt_len = rt_len_arr[rt_len_arr.length>>1];
		// filter by count
		if (rt_len >= opt.min_rt_len) {
			if (s.length < opt.min_cnt_rt) return;
		} else {
			if (s.length < opt.min_cnt) return;
		}
		// count
		let mapq = 0, cnt = {}, cnt_arr = [];
		for (let i = 0; i < s.length; ++i) {
			mapq += s[i]._mapq;
			if (cnt[s[i].source] == null) cnt[s[i].source] = [0, 0];
			cnt[s[i].source][s[i].strand === "+"? 0 : 1]++;
		}
		mapq = (mapq / s.length).toFixed(0);
		let info = `avg_mapq=${mapq};`;
		for (const src in cnt)
			cnt_arr.push(`${src}:${cnt[src][0]},${cnt[src][1]}`);
		info += `count=${cnt_arr.join("|")};`;
		info += `rt_len=${rt_len};`;
		info += v.info.replace(/(;?)source=[^;\s=]+/, "");

		if (!v.is_bp) {
			print(v.ctg, v.st, v.en, ".", s.length, v.strand, info);
		} else {
			print(v.ctg, v.pos, v.ori, v.ctg2, v.pos2, ".", s.length, v.strand, info);
		}
	}

	let last_ctg = null, sv = [];
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		let v = parse_sv(t);
		while (sv.length) {
			if (sv[0].ctg != v.ctg || v.pos - sv[0].pos_max > opt.win_size)
				write_sv(opt, sv.shift().v);
			else break;
		}
		let cnt_same = [];
		for (let i = 0; i < sv.length; ++i) {
			let c = 0;
			if (sv[i].SVTYPE === v.SVTYPE) {
				for (let j = 0; j < sv[i].v.length; ++j)
					if (same_sv(opt, sv[i].v[j], v))
						++c;
			}
			cnt_same[i] = c;
		}
		let max = 0, max_i = -1;
		for (let i = 0; i < sv.length; ++i)
			if (cnt_same[i] > max)
				max = cnt_same[i], max_i = i;
		if (max > 0 && max_i >= 0) { // add to an existing variant
			sv[max_i].v.push(v);
			sv[max_i].pos_max = sv[max_i].pos_max > v.pos? sv[max_i].pos_max : v.pos;
		} else { // create a new variant
			sv.push({ ctg:v.ctg, pos_max:v.pos, SVTYPE:v.SVTYPE, is_bp:v.is_bp, v:[v] });
		}
	}
	while (sv.length)
		write_sv(opt, sv.shift().v);
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: mgutils-es6.js <command> [arguments]");
		print("Commands:");
		print("  merge2vcf    convert merge BED output to VCF");
		print("  addsample    add sample names to merged BED (as a fix)");
		print("  getsv        extract long INDELs and breakpoints from GAF");
		print("  mergesv      merge INDELs and breakpoints from getsv");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'merge2vcf') mg_cmd_merge2vcf(args);
	else if (cmd == 'addsample') mg_cmd_addsample(args);
	else if (cmd == 'getsv') mg_cmd_getsv(args);
	else if (cmd == 'mergesv') mg_cmd_mergesv(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
