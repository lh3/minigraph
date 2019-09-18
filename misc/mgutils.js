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
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'renamefa') mg_cmd_renamefa(args);
	else if (cmd = 'subgaf') mg_cmd_subgaf(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
