#!/usr/bin/env k8

const version = "r605";

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
				for (let i = 5; i < t.length; ++i)
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

function mg_cmd_getlcr(args) {
	let ext = 5, min_lcr = 0.7, min_ac = 5, ref_idx = 0;
	for (const o of getopt(args, "a:r:e:f:")) {
		if (o.opt == "-r") ref_idx = parseInt(o.arg);
		else if (o.opt == "-e") ext = parseInt(o.arg);
		else if (o.opt == "-f") min_lcr = parseFloat(o.arg);
		else if (o.opt == "-a") min_ac = parseInt(o.arg);
	}
	if (args.length == 0) {
		print("Usage: mgutils-es6.js getlcr [options] <merged.bed>");
		print("Options:");
		print(`  -r INT      index of the reference sample [${ref_idx}]`);
		print(`  -f FLOAT    min fraction LCR [${min_lcr}]`);
		print(`  -a INT      min allele count [${min_ac}]`);
		print(`  -e INT      extend by INT-bp on each side in output [${ext}]`);
		return 1;
	}
	const re = /([^\s=;]+)=([^\s=;]+)/g;
	for (const line of k8_readline(args[0])) {
		if (line[0] == '#') continue;
		let m, t = line.split("\t", 5 + ref_idx);
		let ldust = 0, lbb = 0, anno = null, alen = null, ac = null, ref = null;
		while ((m = re.exec(t[3])) != null) {
			if (m[1] == "LBUBBLE") lbb = parseInt(m[2]);
			else if (m[1] == "LDUST") ldust = parseInt(m[2]);
			else if (m[1] == "ANNO") anno = m[2];
			else if (m[1] == "ALEN") alen = m[2].split(",");
			else if (m[1] == "AC") ac = m[2].split(",");
		}
		if (alen == null) continue;
		let is_lcr = /^(lcr|mini|micro|ldust)$/.test(anno);
		if (anno == "segdup" && lbb > 0 && ldust >= lbb * min_lcr)
			is_lcr = true;
		if (!is_lcr) continue;

		if ((m = /^(\d+)/.exec(t[4 + ref_idx])) == null) continue;
		ref = parseInt(m[1]);
		let alen_sel = [];
		for (let i = 0; i < ac.length; ++i) {
			ac[i] = parseInt(ac[i]);
			alen[i] = parseInt(alen[i]);
			if (i == ref || ac[i] >= min_ac)
				alen_sel.push(alen[i]);
		}
		if (alen_sel.length < 2) continue;
		const ctg = t[0].replace(/^[^\s#]+#\d#/, "");
		let st = parseInt(t[1]);
		let en = parseInt(t[2]);
		let max = en - st;
		for (let i = 0; i < alen_sel.length; ++i) {
			const l = parseInt(alen_sel[i]);
			max = l > max? l : max;
		}
		st = st > ext? st - ext : 0;
		print(ctg, st, en + ext, "mg", max);
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
		print("  getlcr       extract LCRs from merged BED");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd === "merge2vcf") mg_cmd_merge2vcf(args);
	else if (cmd === "addsample") mg_cmd_addsample(args);
	else if (cmd === "getlcr") mg_cmd_getlcr(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
