#!/usr/bin/env k8

const gc_version = "r112";

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

// interal query
function iit_sort_copy(a) {
	a.sort((x, y) => (x.st - y.st));
	const b = [];
	for (let i = 0; i < a.length; ++i)
		b.push({ st: a[i].st, en: a[i].en, max: 0, data: a[i].data });
	return b;
}

function iit_index(a) {
	if (a.length == 0) return -1;
	let last, last_i, k;
	for (let i = 0; i < a.length; i += 2) last = a[i].max = a[i].en, last_i = i;
	for (k = 1; 1<<k <= a.length; ++k) {
		const i0 = (1<<k) - 1, step = 1<<(k+1), x = 1<<(k-1);
		for (let i = i0; i < a.length; i += step) {
			a[i].max = a[i].en;
			if (a[i].max < a[i-x].max) a[i].max = a[i-x].max;
			const e = i + x < a.length? a[i+x].max : last;
			if (a[i].max < e) a[i].max = e;
		}
		last_i = last_i>>k&1? last_i - x : last_i + x;
		if (last_i < a.length) last = last > a[last_i].max? last : a[last_i].max;
	}
	return k - 1;
}

function iit_overlap(a, st, en) {
	let h = 0;
	const stack = [], b = [];
	for (h = 0; 1<<h <= a.length; ++h);
	--h;
	stack.push([(1<<h) - 1, h, 0]);
	while (stack.length) {
		const t = stack.pop();
		const x = t[0], h = t[1], w = t[2];
		if (h <= 3) {
			const i0 = x >> h << h;
			let i1 = i0 + (1<<(h+1)) - 1;
			if (i1 >= a.length) i1 = a.length;
			for (let i = i0; i < i1 && a[i].st < en; ++i)
				if (st < a[i].en) b.push(a[i]);
		} else if (w == 0) { // if left child not processed
			stack.push([x, h, 1]);
			const y = x - (1<<(h-1));
			if (y >= a.length || a[y].max > st)
				stack.push([y, h - 1, 0]);
		} else if (x < a.length && a[x].st < en) {
			if (st < a[x].en) b.push(a[x]);
			stack.push([x + (1<<(h-1)), h - 1, 0]);
		}
	}
	return b;
}

function parseNum(s) {
	var m, x = null;
	if ((m = /^(\d*\.?\d*)([mMgGkK]?)/.exec(s)) != null) {
		x = parseFloat(m[1]);
		if (m[2] == 'k' || m[2] == 'K') x *= 1000;
		else if (m[2] == 'm' || m[2] == 'M') x *= 1000000;
		else if (m[2] == 'g' || m[2] == 'G') x *= 1000000000;
	}
	return Math.floor(x + .499);
}

/********************************
 * Extract SVs from GAF/PAF/SAM *
 ********************************/

function mg_revcomp(s) {
	function complement(x) { return { a:'t', t:'a', g:'c', c:'g' }[x] }
	return s.split('').reverse().map(complement).join('');
}

function gc_cmd_extract(args) {
	let opt = { min_mapq:5, min_mapq_end:30, min_frac:0.7, min_len:100, min_aln_len_end:2000, min_aln_len_mid:50, max_cnt_10k:3,
				dbg:false, polyA_pen:5, polyA_drop:100, name:"foo", cen:{} };
	for (const o of getopt(args, "q:Q:l:dc:a:e:m:n:b:", [])) {
		if (o.opt == "-q") opt.min_mapq = parseInt(o.arg);
		else if (o.opt == "-Q") opt.min_mapq_end = parseInt(o.arg);
		else if (o.opt == "-l") opt.min_len = parseInt(o.arg);
		else if (o.opt == "-d") opt.dbg = true;
		else if (o.opt == "-f") opt.min_frac = parseFloat(o.arg);
		else if (o.opt == "-c") opt.max_cnt_10k = parseInt(o.arg);
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
		print("Usage: gafcall.js extract [options] <stable.gaf>");
		print("Options:");
		print(`  -n STR     sample name [${opt.name}]`);
		print(`  -l INT     min INDEL len [${opt.min_len}]`);
		print(`  -f FLOAT   min mapped query fraction [${opt.min_frac}]`);
		print(`  -c INT     max number of long INDELs per 10kb [${opt.max_cnt_10k}]`);
		print(`  -q INT     min mapq [${opt.min_mapq}]`);
		print(`  -Q INT     min mapq for alignment ends [${opt.min_mapq_end}]`);
		print(`  -e INT     min alignment length at ends [${opt.min_aln_len_end}]`);
		print(`  -m INT     min alignment length in the middle [${opt.min_aln_len_mid}]`);
		print(`  -a INT     penalty for non-polyA bases [${opt.polyA_pen}]`);
		print(`  -b FILE    BED for centromeres []`);
		return;
	}

	let re = /(\d+)([=XIDMSHN])/g; // regex for CIGAR
	let re_path = /([><])([^><:\s]+):(\d+)-(\d+)/g; // regex for path/ctg
	let re_ds = /([\+\-\*:])([A-Za-z\[\]0-9]+)/g; // regex for the ds tag
	let re_tsd = /(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?/; // regex for parsing TSD
	let global_qname = "N/A";

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

	function path2ctg(seg, path_off, is_end) {
		let b = [];
		for (let i = 0, k = 0; i < path_off.length; ++i) {
			if (is_end) {
				while (k < seg.length && seg[k].path_en < path_off[i]) ++k;
			} else {
				while (k < seg.length && seg[k].path_en <= path_off[i]) ++k;
			}
			if (k == seg.length) throw Error(`failed to convert path offset to contig offset for read ${global_qname}`);
			const l = path_off[i] - seg[k].path_st;
			if (seg[k].strand > 0)
				b.push({ seg:k, pos:seg[k].ctg_st + l });
			else
				b.push({ seg:k, pos:seg[k].ctg_en - l });
		}
		return b;
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
						a.push({ st:x, en:x,     len:len,  indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:"." });
					else if (op === "D")
						a.push({ st:x, en:x+len, len:-len, indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:"." });
				}
				if (op == "M" || op == "=" || op == "X" || op == "D")
					x += len;
			}
			if (a.length == 0 || a.length > y.qlen * 1e-4 * opt.max_cnt_10k) continue;
			// set stl/enl and str/enr
			for (let i = 0; i < a.length; ++i) {
				a[i].stl = a[i].str = a[i].st;
				a[i].enl = a[i].enr = a[i].en;
			}
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
							if (a[i].st != x || a[i].en != x || a[i].len != len)
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
					const llen = m[2]? m[2].length : 0;
					const rlen = m[5]? m[5].length : 0;
					a[i].stl = a[i].st - rlen;
					a[i].enl = a[i].en - rlen;
					a[i].str = a[i].st + llen;
					a[i].enr = a[i].en + llen;
				}
			} // ~if(y.ds)
			if (opt.dbg) print('X0', line);
			let seg = []; // reference segments in the path
			if (/[><]/.test(y.path)) { // with ><: this is a path
				let x = 0;
				if (y.strand != '+') throw Error("reverse strand on path");
				while ((m = re_path.exec(y.path)) != null) {
					const st = parseInt(m[3]), en = parseInt(m[4]);
					const strand = m[1] === '>'? 1 : -1;
					seg.push({ ctg:m[2], ctg_st:st, ctg_en:en, strand:strand, path_st:x, path_en:x + (en - st) });
					x += en - st;
				}
			} else { // this is a contig name
				seg.push({ ctg:y.path, ctg_st:0, ctg_en:y.tlen, strand:1, path_st:0, path_en:y.tlen });
			}
			let off_stl = [], off_str = [], off_enl = [], off_enr = [];
			for (let i = 0; i < a.length; ++i) {
				off_stl[i] = a[i].stl, off_enl[i] = a[i].enl;
				off_str[i] = a[i].str, off_enr[i] = a[i].enr;
			}
			global_qname = y.qname;
			const stl = path2ctg(seg, off_stl, false);
			const enl = path2ctg(seg, off_enl, true);
			const str = path2ctg(seg, off_str, false);
			const enr = path2ctg(seg, off_enr, true);
			for (let i = 0; i < a.length; ++i) {
				if (!(stl[i].seg == str[i].seg && stl[i].seg == enl[i].seg && str[i].seg == enr[i].seg)) continue; // all on the same segment
				const s = seg[stl[i].seg];
				let st = stl[i].pos, en = enl[i].pos, strand = y.strand;
				if (s.strand < 0) { // then reverse complement tsd, polyA and insert
					a[i].polyA_len = -a[i].polyA_len;
					a[i].tsd_seq = mg_revcomp(a[i].tsd_seq);
					a[i].int_seq = mg_revcomp(a[i].int_seq);
					st = enr[i].pos, en = str[i].pos;
					strand = y.strand === "+"? "-" : "+";
				}
				let info1 = (a[i].len > 0? "SVTYPE=INS" : "SVTYPE=DEL") + `;SVLEN=${a[i].len};tsd_len=${a[i].tsd_len};polyA_len=${a[i].polyA_len}`;
				const info2 = `source=${opt.name};tsd_seq=${a[i].tsd_seq.length>0?a[i].tsd_seq:"."};insert=${a[i].int_seq.length>0?a[i].int_seq:"."}`;
				if (opt.cen[s.ctg] != null) {
					const dist_st = cal_cen_dist(opt, s.ctg, st);
					const dist_en = cal_cen_dist(opt, s.ctg, en);
					info1 += `;cen_dist=${dist_st < dist_en? dist_st : dist_en}`
				}
				print(s.ctg, st, en, y.qname, y.mapq, strand, `${info1};${info2}`);
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

/*************************
 * Merge extracted calls *
 *************************/

function gc_cmd_merge(args) {
	let opt = { min_cnt:4, min_cnt_strand:2, min_cnt_rt:1, min_rt_len:0, win_size:100, max_diff:0.05, min_cen_dist:500000, max_allele:100, max_check:500 };
	for (const o of getopt(args, "w:d:c:e:r:R:s:A:C:")) {
		if (o.opt === "-w") opt.win_size = parseInt(o.arg);
		else if (o.opt === "-d") opt.max_diff = parseFloat(o.arg);
		else if (o.opt === "-c") opt.min_cnt = parseInt(o.arg);
		else if (o.opt === "-s") opt.min_cnt_strand = parseInt(o.arg);
		else if (o.opt === "-e") opt.min_cen_dist = parseNum(o.arg);
		else if (o.opt === "-r") opt.min_rt_len = parseInt(o.arg);
		else if (o.opt === "-R") opt.min_cnt_rt = parseInt(o.arg);
		else if (o.opt === "-A") opt.max_allele = parseInt(o.arg);
		else if (o.opt === "-C") opt.max_check = parseInt(o.arg);
	}
	if (args.length == 0) {
		print("Usage: sort -k1,1 -k2,2n extract.output | gafcall.js merge [options] -");
		print("Options:");
		print(`  -c INT     min read count [${opt.min_cnt}]`);
		print(`  -s INT     min read count on each strand [${opt.min_cnt_strand}]`);
		print(`  -w INT     window size [${opt.win_size}]`);
		print(`  -d FLOAT   max allele length difference ratio [${opt.max_diff}]`);
		print(`  -e NUM     min distance to centromeres [${opt.min_cen_dist}]`);
		print(`  -r INT     min min(TSD,polyA) to tag a candidate RT; 0 to disable [${opt.min_rt_len}]`);
		print(`  -R INT     min read count for a candidate RT [${opt.min_cnt_rt}]`);
		print(`  -A INT     check up to INT nearby alleles [${opt.max_allele}]`);
		print(`  -C INT     compare up to INT reads per allele [${opt.max_check}]`);
		return;
	}

	function splitmix32(a) { // random number generator as we can't set seeds in Javascript
		return function() {
			a |= 0; a = a + 0x9e3779b9 | 0;
			let t = a ^ a >>> 16;
			t = Math.imul(t, 0x21f0aaad);
			t = t ^ t >>> 15;
			t = Math.imul(t, 0x735a2d97);
			return ((t = t ^ t >>> 15) >>> 0) / 4294967296;
		}
	}

	function parse_sv(t) {
		const re_info = /([^;\s=]+)=([^;\s=]+)/g;
		let v = { ctg:t[0], pos:parseInt(t[1]), st:-1, en:-1, ori:null, ctg2:null, pos2:null, _mapq:0, strand:null, is_bp:false, name:null, info:null }
		v.is_bp = /[><]/.test(t[2]);
		const off = v.is_bp? 6 : 4; // offset of the mapq column
		v._mapq = parseInt(t[off]);
		v.strand = t[off+1];
		v.info = t[off+2];
		v.name = t[off-1];
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
			if (!((v.ori === "><" && w.ori === "<>") || (v.ori === "<>" && w.ori === "><")))
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
		// count
		let mapq = 0, cnt = {}, cnt_strand = [0, 0], name = [];
		for (let i = 0; i < s.length; ++i) {
			mapq += s[i]._mapq;
			if (cnt[s[i].source] == null) cnt[s[i].source] = [0, 0];
			cnt[s[i].source][s[i].strand === "+"? 0 : 1]++;
			cnt_strand[s[i].strand === "+"? 0 : 1]++;
			name.push(s[i].name);
		}
		mapq = (mapq / s.length).toFixed(0);
		// filter by count
		if (opt.min_rt_len > 0 && rt_len >= opt.min_rt_len) {
			if (s.length < opt.min_cnt_rt) return;
		} else {
			if (s.length < opt.min_cnt) return;
			if (cnt_strand[0] < opt.min_cnt_strand || cnt_strand[1] < opt.min_cnt_strand) return;
		}
		let info = `avg_mapq=${mapq};`, cnt_arr = [];
		for (const src in cnt)
			cnt_arr.push(`${src}:${cnt[src][0]},${cnt[src][1]}`);
		info += `count=${cnt_arr.join("|")};`;
		info += `rt_len=${rt_len};`;
		info += v.info.replace(/(;?)source=[^;\s=]+/, "");
		info += `;reads=${name.join(",")}`;

		if (!v.is_bp) {
			print(v.ctg, v.st, v.en, ".", s.length, v.strand, info);
		} else {
			print(v.ctg, v.pos, v.ori, v.ctg2, v.pos2, ".", s.length, v.strand, info);
		}
	}

	let last_ctg = null, sv = [], rng = splitmix32(11);
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		let v = parse_sv(t);
		while (sv.length) {
			if (sv[0].ctg != v.ctg || v.pos - sv[0].pos_max > opt.win_size || sv.length > opt.max_allele)
				write_sv(opt, sv.shift().v);
			else break;
		}
		let cnt_same = [];
		for (let i = 0; i < sv.length; ++i) {
			let c = 0;
			if (sv[i].SVTYPE === v.SVTYPE || (sv[i].SVTYPE === "INS" && v.SVTYPE === "DUP")) {
				if (sv[i].v.length <= opt.max_check) {
					for (let j = 0; j < sv[i].v.length; ++j)
						if (same_sv(opt, sv[i].v[j], v))
							++c;
				} else { // use reservior sampling to samplg a subset of reads
					let p = [];
					for (let j = 0; j < sv[i].v.length; ++j) {
						let k = j < opt.max_check? j : Math.floor(j * rng());
						if (k < opt.max_check) p[k] = j;
					}
					for (let k = 0; k < opt.max_check; ++k)
						if (same_sv(opt, sv[i].v[p[k]], v))
							++c;
					c = Math.floor(c / opt.max_check * sv[i].v.length + .499);
				}
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
			sv.push({ ctg:v.ctg, pos_max:v.pos, SVTYPE: v.SVTYPE != "DUP"? v.SVTYPE : "INS", is_bp:v.is_bp, v:[v] });
		}
	}
	while (sv.length)
		write_sv(opt, sv.shift().v);
}

/*************************
 * Parse and reformat SV *
 *************************/

function gc_parse_sv(min_len, fn, ignore_flt, check_gt) {
	let sv = [], ignore_id = {};
	ignore_flt = typeof ignore_flt !== "undefined"? ignore_flt : true;
	check_gt = typeof check_gt !== "undefined"? check_gt : false;
	for (const line of k8_readline(fn)) {
		if (line[0] === "#") continue;
		let m, t = line.split("\t");
		if (!/^\d+$/.test(t[1])) continue;
		//print("X", line);
		t[1] = parseInt(t[1]);
		let type = 0, info = null;
		if (/^[><][><]$/.test(t[2])) type = 3, info = t[8]; // breakpoint
		else if (/;/.test(t[7])) type = 1, info = t[7]; // VCF
		else if (/^\d+$/.test(t[2]) && /;/.test(t[6])) type = 2, info = t[6]; // BED
		if (type == 0) continue;
		let svtype = null, svlen = 0;
		if ((m = /SVTYPE=([^\s;]+)/.exec(info)) != null)
			svtype = m[1];
		if ((m = /SVLEN=([^\s;]+)/.exec(info)) != null)
			svlen = parseFloat(m[1]);
		if (type == 2) { // BED line
			t[2] = parseInt(t[2]);
			if (t[1] > t[2]) throw("incorrected BED?");
			if (Math.abs(svlen) < min_len) continue;
			sv.push({ ctg:t[0], pos:t[1], ctg2:t[0], pos2:t[2], ori:">>", svtype:svtype, svlen:svlen, vaf:1 });
		} else if (type == 3) { // breakpoint line
			t[4] = parseInt(t[4]);
			if (t[0] === t[3] && Math.abs(svlen) < min_len) continue;
			sv.push({ ctg:t[0], pos:t[1], ctg2:t[3], pos2:t[4], ori:t[2], svtype:svtype, svlen:svlen, vaf:1 });
		} else if (type == 1) { // VCF line
			if (!ignore_flt && t[6] !== "PASS" && t[6] !== ".") continue; // ignore filtered calls
			if (check_gt && t.length >= 9 && /^0[\/\|]0/.test(t[9])) continue; // not a variant
			let rlen = t[3].length, en = t[1] + rlen - 1;
			let s = { ctg:t[0], pos:t[1]-1, ctg2:t[0], pos2:en, ori:">>", vaf:1 };
			if ((m = /\bVAF=([^\s;]+)/.exec(info)) != null)
				s.vaf = parseFloat(m[1]);
			if (/^[A-Z,\*]+$/.test(t[4])) { // assume full allele sequence; override SVTYPE/SVLEN even if present
				let alt = t[4].split(",");
				for (let i = 0; i < alt.length; ++i) {
					const a = alt[i], len = a.length - rlen;
					if (Math.abs(len) < min_len) continue;
					if (len < 0)
						sv.push({ ctg:s.ctg, pos:s.pos, ctg2:s.ctg, pos2:en, svtype:"DEL", svlen:len, ori:">>", vaf:s.vaf });
					else
						sv.push({ ctg:s.ctg, pos:s.pos, ctg2:s.ctg, pos2:en, svtype:"INS", svlen:len, ori:">>", vaf:s.vaf });
				}
			} else { // other SV encoding
				if (t[2] !== ".") {
					if (ignore_id[t[2]]) continue; // ignore previously visited ID
					ignore_id[t[2]] = 1;
				}
				if ((m = /\b(MATE_ID|MATEID)=([^\s;]+)/.exec(info)) != null)
					ignore_id[m[2]] = 1;
				if (svtype == null) throw Error(`can't determine SVTYPE: ${t.join("\t")}`); // we don't infer SVTYPE from breakpoint
				s.svtype = svtype;
				if (svtype !== "BND" && Math.abs(svlen) < min_len) continue; // too short
				if (svtype === "DEL" && svlen > 0) svlen = -svlen; // correct SVLEN as some VCF encodes this differently
				s.svlen = svlen;
				if ((m = /\bEND=(\d+)/.exec(info)) != null) {
					s.pos2 = parseInt(m[1]);
				} else if (rlen == 1) {
					if (svtype === "BND" && t[4].length < 6) continue; // ignore one-sided breakpoint
					if (svtype === "DEL" || svtype === "DUP" || svtype === "INV")
						s.pos2 = s.pos + Math.abs(svlen);
				}
				if ((m = /^[A-Z]+\[([^\s:]+):(\d+)\[$/.exec(t[4])) != null) s.ctg2 = m[1], s.pos2 = parseInt(m[2]), s.ori = ">>";
				else if ((m = /^\]([^\s:]+):(\d+)\][A-Z]+$/.exec(t[4])) != null) s.ctg2 = m[1], s.pos2 = parseInt(m[2]), s.ori = "<<";
				else if ((m = /^\[([^\s:]+):(\d+)\[[A-Z]+$/.exec(t[4])) != null) s.ctg2 = m[1], s.pos2 = parseInt(m[2]), s.ori = "<>";
				else if ((m = /^[A-Z]+\]([^\s:]+):(\d+)\]$/.exec(t[4])) != null) s.ctg2 = m[1], s.pos2 = parseInt(m[2]), s.ori = "><";
				if (svtype !== "BND" && s.ctg !== s.ctg2) throw Error("different contigs for non-BND type");
				if (svtype === "BND" && s.ctg === s.ctg2 && Math.abs(svlen) < min_len) continue;
				if (s.ctg === s.ctg2 && s.pos > s.pos2) {
					let tmp = s.pos;
					s.pos = s.pos2, s.pos2 = tmp;
				}
				sv.push(s);
			}
		}
	}
	return sv;
}

function gc_read_bed(fn) {
	let h = {};
	for (const line of k8_readline(fn)) {
		let t = line.split("\t");
		if (t.length < 3) continue;
		if (h[t[0]] == null) h[t[0]] = [];
		h[t[0]].push({ st:parseInt(t[1]), en:parseInt(t[2]), data:null });
	}
	for (const ctg in h) {
		h[ctg] = iit_sort_copy(h[ctg]);
		iit_index(h[ctg]);
	}
	return h;
}

function gc_cmd_view(args) {
	let min_read_len = 100, ignore_flt = false, check_gt = false, count_long = false, bed = null;
	for (const o of getopt(args, "l:FGCb:")) {
		if (o.opt === "-l") min_read_len = parseNum(o.arg);
		else if (o.opt === "-F") ignore_flt = true;
		else if (o.opt === "-G") check_gt = true;
		else if (o.opt === "-C") count_long = true;
		else if (o.opt === "-b") bed = gc_read_bed(o.arg);
	}
	if (args.length == 0) {
		print("Usage: gafcall.js view [options] <in.vcf>");
		print("Options:");
		print(`  -l NUM       min length [${min_read_len}]`);
		print(`  -b FILE      regions to include []`);
		print(`  -F           ignore FILTER field in VCF`);
		print(`  -G           check GT in VCF`);
		print(`  -C           count 20kb, 100kb, 1Mb and translocations`);
		return;
	}
	for (let j = 0; j < args.length; ++j) {
		const sv = gc_parse_sv(min_read_len, args[j], ignore_flt, check_gt);
		let cnt = [ 0, 0, 0, 0 ];
		for (let i = 0; i < sv.length; ++i) {
			const s = sv[i];
			if (bed != null) {
				if (bed[s.ctg] == null || bed[s.ctg2] == null) continue;
				if (iit_overlap(bed[s.ctg],  s.pos,  s.pos  + 1).length === 0) continue;
				if (iit_overlap(bed[s.ctg2], s.pos2, s.pos2 + 1).length === 0) continue;
			}
			if (count_long) {
				if (s.ctg != s.ctg2) {
					++cnt[0], ++cnt[1], ++cnt[2], ++cnt[3];
				} else {
					const len = Math.abs(s.svlen);
					if (len >= 1000000) ++cnt[1];
					if (len >= 100000) ++cnt[2];
					if (len >= 20000) ++cnt[3];
				}
			} else {
				print(s.ctg, s.pos, s.ori, s.ctg2, s.pos2, s.svtype, s.svlen);
			}
		}
		if (count_long)
			print(cnt.join("\t"), args[j]);
	}
}

/**************
 * Evaluation *
 **************/

function gc_cmp_sv(opt, base, test, label) {
	let h = {};
	for (let i = 0; i < base.length; ++i) {
		const s = base[i];
		if (h[s.ctg] == null) h[s.ctg] = [];
		if (h[s.ctg2] == null) h[s.ctg2] = [];
		h[s.ctg].push({ st:s.pos, en:s.pos+1, data:s });
		h[s.ctg2].push({ st:s.pos2, en:s.pos2+1, data:s });
	}
	for (const ctg in h) {
		h[ctg] = iit_sort_copy(h[ctg]);
		iit_index(h[ctg]);
	}

	function same_sv1(opt, b, t) { // compare two SVs
		// check type
		if (b.svtype != t.svtype) { // type mismatch
			if (!(b.svtype === "DUP" && t.svtype === "INS") && !(b.svtype === "INS" && t.svtype === "DUP") && b.svtype !== "BND" && t.svtype !== "BND") // special case for INS vs DUP
				return false;
		}
		// check length
		const len_check = (Math.abs(b.svlen) >= Math.abs(t.svlen) * opt.min_len_ratio && Math.abs(t.svlen) >= Math.abs(b.svlen) * opt.min_len_ratio);
		if (b.svtype !== "BND" && t.svtype !== "BND" && !len_check) return false;
		// check the coordinates of end points
		let match1 = 0, match2 = 0;
		if (t.ctg == b.ctg   && t.pos >= b.pos - opt.win_size   && t.pos <= b.pos + opt.win_size)   match1 |= 1;
		if (t.ctg == b.ctg2  && t.pos >= b.pos2 - opt.win_size  && t.pos <= b.pos2 + opt.win_size)  match1 |= 2;
		if (t.ctg2 == b.ctg  && t.pos2 >= b.pos - opt.win_size  && t.pos2 <= b.pos + opt.win_size)  match2 |= 1;
		if (t.ctg2 == b.ctg2 && t.pos2 >= b.pos2 - opt.win_size && t.pos2 <= b.pos2 + opt.win_size) match2 |= 2;
		if (b.svtype === "DUP" && t.svtype === "INS") {
			return ((match1&1) != 0);
		} else if (b.svtype === "INS" && t.svtype === "DUP") {
			return ((match1&1) != 0);
		} else if (b.svtype === "BND" || t.svtype === "BND") {
			return (((match1&1) != 0 && (match2&2) != 0) || ((match1&2) != 0 && (match2&1) != 0));
		} else {
			return ((match1&1) != 0 && (match2&2) != 0);
		}
	}

	function eval1(opt, h, ctg, pos, t) {
		if (h[ctg] == null) return false;
		const st = pos > opt.win_size? pos - opt.win_size : 0;
		const en = pos + opt.win_size;
		const a = iit_overlap(h[ctg], st, en);
		let n = 0;
		for (let i = 0; i < a.length; ++i)
			if (same_sv1(opt, a[i].data, t))
				++n;
		return n;
	}

	let tot = 0, error = 0;
	for (let j = 0; j < test.length; ++j) {
		const t = test[j];
		if (t.svtype !== "BND" && Math.abs(t.svlen) < opt.min_len) continue; // not long enough for non-BND type; note that t.ctg === t.ctg2 MUST stand due to assertion in parsing
		if (t.svtype === "BND" && t.ctg === t.ctg2 && Math.abs(t.svlen) < opt.min_len) continue; // not long enough; in principle, this can be merged to the line above
		if (t.vaf != null && t.vaf < opt.min_vaf) continue; // filter by VAF
		if (opt.bed != null) {
			if (opt.bed[t.ctg] == null || opt.bed[t.ctg2] == null) continue;
			if (iit_overlap(opt.bed[t.ctg],  t.pos,  t.pos  + 1).length === 0) continue;
			if (iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2 + 1).length === 0) continue;
		}
		++tot;
		const n = eval1(opt, h, t.ctg, t.pos, t) + eval1(opt, h, t.ctg2, t.pos2, t);
		if (n == 0) {
			++error;
			if (opt.print_err)
				print(label, t.ctg, t.pos, t.ori, t.ctg2, t.pos2, t.svtype, t.svlen);
		}
	}
	return [tot, error];
}

function gc_cmd_eval(args) {
	let opt = { min_len:100, read_len_ratio:0.8, win_size:500, min_len_ratio:0.6, min_vaf:0, bed:null, dbg:false, print_err:false, ignore_flt:false, check_gt:false };
	for (const o of getopt(args, "dr:l:w:em:v:b:FG")) {
		if (o.opt === "-d") opt.dbg = true;
		else if (o.opt === "-b") opt.bed = gc_read_bed(o.arg);
		else if (o.opt === "-l") opt.min_len = parseNum(o.arg);
		else if (o.opt === "-m") opt.min_len_ratio = parseFloat(o.arg);
		else if (o.opt === "-r") opt.read_len_ratio = parseFloat(o.arg);
		else if (o.opt === "-w") opt.win_size = parseNum(o.arg);
		else if (o.opt === "-v") opt.min_vaf = parseFloat(o.arg);
		else if (o.opt === "-F") opt.ignore_flt = true;
		else if (o.opt === "-G") opt.check_gt = true;
		else if (o.opt === "-e") opt.print_err = true;
	}
	if (args.length < 2) {
		print("Usgae: gafcall.js eval [options] <base.vcf> <test.vcf>");
		print("Options:");
		print(`  -b FILE     confident regions in BED []`);
		print(`  -l NUM      min SVLEN [${opt.min_len}]`);
		print(`  -w NUM      fuzzy window size [${opt.win_size}]`);
		print(`  -r FLOAT    read SVs longer than {-l}*FLOAT [${opt.read_len_ratio}]`);
		print(`  -m FLOAT    two SVs regarded the same if length ratio above [${opt.min_len_ratio}]`);
		print(`  -v FLOAT    ignore VAF below FLOAT (requiring VAF in VCF) [${opt.min_vaf}]`);
		print(`  -F          ignore FILTER in VCF`);
		print(`  -G          check GT in VCF`);
		print(`  -e          print errors`);
		return;
	}
	const min_read_len = Math.floor(opt.min_len * opt.read_len_ratio + .499);

	if (args.length === 2) { // two-sample mode
		const base = gc_parse_sv(min_read_len, args[0], opt.ignore_flt, opt.check_gt);
		const test = gc_parse_sv(min_read_len, args[1], opt.ignore_flt, opt.check_gt);
		const [tot_fn, fn] = gc_cmp_sv(opt, test, base, "FN");
		const [tot_fp, fp] = gc_cmp_sv(opt, base, test, "FP");
		print("RN", tot_fn, fn, (fn / tot_fn).toFixed(4), args[0]);
		print("RP", tot_fp, fp, (fp / tot_fp).toFixed(4), args[1]);
	} else { // multi-sample mode
		let vcf = [];
		for (let i = 0; i < args.length; ++i)
			vcf[i] = gc_parse_sv(min_read_len, args[i], opt.ignore_flt, opt.check_gt);
		for (let i = 0; i < args.length; ++i) {
			let a = [ "SN" ];
			for (let j = 0; j < args.length; ++j) {
				const [cnt, err] = gc_cmp_sv(opt, vcf[i], vcf[j], "XX");
				if (i != j) a.push((1 - err/cnt).toFixed(4));
				else a.push(cnt);
			}
			print(a.join("\t"), args[i]);
		}
	}
}

/**********************
 * Join two GSV files *
 **********************/

function gc_cmd_join(args) {
	for (const o of getopt(args, "")) {
	}
	if (args.length < 2) {
		print("Usgae: gafcall.js join <filter.gsv> <out.gsv>");
		return;
	}

	function get_type(t, col_info) {
		const info = t[col_info];
		let m;
		if ((m = /\bSVTYPE=([^\s;]+)/.exec(info)) != null) {
			if (m[1] === "INS" || m[1] === "DUP") return 1;
			else if (m[1] === "DEL") return 2;
			else if (m[1] === "INV") return 4;
			else if (m[1] === "BND" && col_info === 8 && t[0] !== t[3]) return 8;
		}
		return 0;
	}

	let h = {}
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		const col_info = /^[><]+$/.test(t[2])? 8 : 6;
		const name = t[col_info - 3];
		if (h[name] == null) h[name] = 0;
		h[name] |= get_type(t, col_info);
	}
	for (const line of k8_readline(args[1])) {
		let t = line.split("\t");
		const col_info = /^[><]+$/.test(t[2])? 8 : 6;
		const name = t[col_info - 3];
		if (h[name] == null) continue;
		const type = get_type(t, col_info);
		if (type === 0 || (h[name]&type) || (h[name]&8))
			print(line);
	}
}

/*******************************
 * Convert to VCF (unfinished) *
 *******************************/

function mg_cmd_sv2vcf(args) {
	let opt = { };
	for (const o of getopt(args, "")) {
	}
	if (args.length == 0) {
		print("Usage: mgutils-es6.js sv2vcf [options] <sv>");
		return;
	}

	const re_info = /([^;\s=]+)=([^;\s=]+)/g;
	const key = { SVTYPE:1, SVLEN:1 };
	print(`##fileformat=VCFv4.3`);
	print(`##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">`);
	print(`##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">`);
	print(`##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">`);
	print(`##ALT=<ID=DEL,Description="Deletion">`);
	print(`##ALT=<ID=INS,Description="Insertion">`);
	print(`##ALT=<ID=DUP,Description="Duplication">`);
	print(`##ALT=<ID=INV,Description="Inversion">`);
	print(`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`);
	print(`#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample`);
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		const is_bp = /[><]/.test(t[2]);
		if (is_bp && t[0] !== t[3]) continue;
		const off_info = is_bp? 8 : 6;
		let m, type = null, info = "";
		while ((m = re_info.exec(t[off_info])) != null) {
			if (key[m[1]]) {
				if (info.length) info += `;`;
				info += `${m[1]}=${m[2]}`;
			}
			if (m[1] === "SVTYPE") type = m[2];
		}
		if (type == null || type === "BND") continue;
		info += is_bp? `;END=${t[4]}` : `;END=${t[2]}`;
		print(t[0], t[1], ".", "N", `<${type}>`, t[off_info-2], `.`, info, "GT", "1/1");
	}
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: gafcall.js <command> [arguments]");
		print("Commands:");
		print("  extract      extract long INDELs and breakpoints from GAF");
		print("  merge        merge extracted INDELs and breakpoints");
		print("  eval         evaluate SV calls");
		print("  view         print in the gafcall format");
		print("  join         join two 'extract' outputs");
		print("  version      print version number");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd === "extract" || cmd === "getsv") gc_cmd_extract(args);
	else if (cmd === "merge" || cmd === "mergesv") gc_cmd_merge(args);
	else if (cmd === "eval") gc_cmd_eval(args);
	else if (cmd === "view" || cmd === "format") gc_cmd_view(args);
	else if (cmd === "join") gc_cmd_join(args);
	else if (cmd === "version") {
		print(gc_version);
		return;
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
