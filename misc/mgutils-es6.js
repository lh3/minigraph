#!/usr/bin/env k8

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
		warn("Usage: mgutils-es6.js merge2vcf [options] <in.bed>");
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
		hdr.push(`##ALT=<ID=CNV:${i},Description="Allele ${i}">`);
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

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: mgutils-es6.js <command> [arguments]");
		print("Commands:");
		print("  merge2vcf    convert merge BED output to VCF");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'merge2vcf') mg_cmd_merge2vcf(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
