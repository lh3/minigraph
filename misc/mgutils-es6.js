#!/usr/bin/env k8

const version = "r565";

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

function mg_cmd_getindel(args) {
	let min_mapq = 5, min_frac = 0.7, min_len = 100, max_cnt = 5, dbg = false;
	for (const o of getopt(args, "q:l:dc:", [])) {
		if (o.opt == "-q") min_mapq = parseInt(o.arg);
		else if (o.opt == "-l") min_len = parseInt(o.arg);
		else if (o.opt == "-d") dbg = true;
		else if (o.opt == "-f") min_frac = parseFloat(o.arg);
		else if (o.opt == "-c") max_cnt = parseInt(o.arg);
	}
	if (args.length == 0) {
		print("Usage: mgutils-es6.js getindel [options] <stable.gaf>");
		print("Options:");
		print(`  -q INT     min mapq [${min_mapq}]`);
		print(`  -l INT     min INDEL len [${min_len}]`);
		print(`  -f FLOAT   min mapped query fraction [${min_frac}]`);
		print(`  -c INT     max count per read [${max_cnt}]`);
		return;
	}
	let re = /(\d+)([=XIDM])/g, re_path = /([><])([^><:\s]+):(\d+)-(\d+)/g;
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		if (t.length < 12) continue;
		for (let i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
		for (let i = 6; i <= 11; ++i) t[i] = parseInt(t[i]);
		if (t[11] < min_mapq) continue;
		if (t[3] - t[2] < t[1] * min_frac) continue;
		let cg = null;
		for (let i = 12; i < t.length; ++i)
			if (t[i].substr(0, 5) === "cg:Z:")
				cg = t[i].substr(5);
		if (cg == null) continue;
		let m, a = [], x = t[7];
		while ((m = re.exec(cg)) != null) {
			const op = m[2], len = parseInt(m[1]);
			if (len >= min_len) {
				if (op === "I") a.push([x - 1, x + 1, len]);
				else if (op === "D") a.push([x, x + len, -len]);
			}
			if (op == "M" || op == "=" || op == "X" || op == "D")
				x += len;
		}
		if (a.length == 0 || a.length > max_cnt) continue;
		if (dbg) print('X0', line);
		let seg = [];
		if (/[><]/.test(t[5])) { // with ><
			let y = 0;
			if (t[4] != '+') throw Error("reverse strand on path");
			while ((m = re_path.exec(t[5])) != null) {
				const st = parseInt(m[3]), en = parseInt(m[4]);
				seg.push([m[2], st, en, m[1] == '>'? 1 : -1, y, y + (en - st)]);
				y += en - st;
			}
		} else {
			seg.push([t[5], 0, t[6], 1, 0, t[6]]);
		}
		let st = [], en = [];
		for (let i = 0, k = 0; i < a.length; ++i) { // start
			while (k < seg.length && seg[k][5] <= a[i][0]) ++k;
			if (k == seg.length) throw Error("failed to find start");
			const l = a[i][0] - seg[k][4];
			if (seg[k][3] > 0)
				st.push([k, seg[k][1] + l, seg[k][1] + l]);
			else
				st.push([k, seg[k][2] - l, seg[k][1] + l]);
		}
		for (let i = 0, k = 0; i < a.length; ++i) { // end
			while (k < seg.length && seg[k][5] <= a[i][1]) ++k;
			if (k == seg.length) throw Error("failed to find end");
			const l = a[i][1] - seg[k][4];
			if (seg[k][3] > 0)
				en.push([k, seg[k][1] + l, seg[k][1] + l]);
			else
				en.push([k, seg[k][2] - l, seg[k][1] + l]);
		}
		for (let i = 0; i < a.length; ++i) {
			if (dbg) print('X2', a[i][0], a[i][1], st[i][0], en[i][0]);
			if (st[i][0] == en[i][0]) { // on the same segment
				const s = seg[st[i][0]];
				const strand = s[3] > 0? t[4] : t[4] == '+'? '-' : '+';
				print(s[0], st[i][1] < en[i][1]? st[i][1] : en[i][1], st[i][1] > en[i][1]? st[i][1] : en[i][1], t[0], t[11], strand, a[i][2]);
			} else { // on different segments
				let path = [], len = 0;
				for (let j = st[i][0]; j <= en[i][0]; ++j) {
					const s = seg[j];
					len += s[2] - s[1];
					path.push((s[3] > 0? '>' : '<') + s[0] + `:${s[1]}-${s[2]}`);
				}
				const off = seg[st[i][0]][4];
				print(path.join(""), a[i][0] - off, a[i][1] - off, t[0], t[11], '+', a[i][2]);
			}
		}
	}
}

function mg_cmd_mergeindel(args) {
	let min_mapq = 5, min_cnt = 1, win_size = 100, max_diff = 0.05;
	for (const o of getopt(args, "q:c:w:d:", [])) {
		if (o.opt == "-q") min_mapq = parseInt(o.arg);
		else if (o.opt == "-c") min_cnt = parseInt(o.arg);
		else if (o.opt == "-w") win_size = parseInt(o.arg);
		else if (o.opt == "-d") max_diff = parseFloat(o.arg);
	}
	if (args.length == 0) {
		print("Usage: mgutils-es6.js mergeindel [options] <stable.gaf>");
		print("Options:");
		print(`  -q INT     min average mapq [${min_mapq}]`);
		print(`  -c INT     min read count [${min_cnt}]`);
		print(`  -w INT     window size [${win_size}]`);
		print(`  -d FLOAT   max allele length difference ratio [${max_diff}]`);
		return;
	}
	let h = {};
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		t[1] = parseInt(t[1]);
		t[2] = parseInt(t[2]);
		t[4] = parseInt(t[4]);
		t[6] = parseInt(t[6]);
		const ctg = t.shift();
		if (h[ctg] == null) h[ctg] = [];
		h[ctg].push(t);
	}

	function print_bed(ctg, t) {
		let len = 0, mapq = 0, n = t[6].length, nf = 0, nr = 0;
		if (n < min_cnt) return;
		for (let i = 0; i < n; ++i) {
			len += t[6][i];
			mapq += t[7][i];
			if (t[8][i][0] == '+') ++nf;
			else ++nr;
		}
		len = Math.floor(len / n + .499);
		const len_str = len > 0? "+" + len.toString() : len.toString();
		mapq = Math.floor(mapq / n + .499);
		if (mapq < min_mapq) return;
		print(ctg, t[0], t[1], len_str, n, ".", `mq:i:${mapq}`, `cf:i:${nf}`, `cr:i:${nr}`, `rd:Z:${t[8].join(",")}`);
	}

	for (const ctg in h) {
		h[ctg].sort(function(x,y) { return x[0]-y[0]; });
		const a = h[ctg];
		let b = [];
		for (let i = 0; i < a.length; ++i) {
			const ai = a[i];
			while (b.length) {
				if (ai[0] - b[0][1] > win_size) {
					const t = b.shift();
					print_bed(ctg, t);
				} else break;
			}
			let merge_j = -1;
			for (let j = b.length - 1; j >= 0; --j) {
				let bj = b[j];
				if (bj[5] * ai[5] <= 0) continue; // not the same sign
				const la = ai[5] > 0? ai[5] : -ai[5];
				const lb = bj[5] > 0? bj[5] : -bj[5];
				const diff = la > lb? la - lb : lb - la;
				if (diff > (la > lb? la : lb) * max_diff) continue;
				bj[6].push(ai[5]); // length
				bj[7].push(ai[3]); // mapQ
				bj[8].push(`${ai[4]}${ai[2]}`); // read name
				bj[1] = bj[1] > ai[1]? bj[1] : ai[1];
				merge_j = j;
				break;
			}
			if (merge_j < 0)
				b.push([ai[0], ai[1], ".", ai[3], ".", ai[5], [ai[5]], [ai[3]], [`${ai[4]}${ai[2]}`]]);
		}
		while (b.length) {
			const t = b.shift();
			print_bed(ctg, t);
		}
	}
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
		print("  getindel     extract long INDELs from GAF");
		print("  mergeindel   merge long INDELs from getindel");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'merge2vcf') mg_cmd_merge2vcf(args);
	else if (cmd == 'addsample') mg_cmd_addsample(args);
	else if (cmd == 'getindel') mg_cmd_getindel(args);
	else if (cmd == 'mergeindel') mg_cmd_mergeindel(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
