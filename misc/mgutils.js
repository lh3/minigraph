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

function mg_cmd_fa2gfa(args)
{
	var c, rank = 0;
	while ((c = getopt(args, "r:")) != null)
		if (c == 'r') rank = parseInt(getopt.arg);
	if (getopt.ind == args.length) {
		print("Usage: mgutils.js fa2gfa [-r rank] <in.fa>");
		return;
	}
	var file = new File(args[getopt.ind]);
	var buf = new Bytes();
	var seq = new Bytes();
	var name = null, id = 1;
	while (file.readline(buf) >= 0) {
		if (buf[0] != 62) {
			seq.set(buf);
		} else {
			var m, s = buf.toString();
			if ((m = /^>(\S+)/.exec(s)) != null) {
				if (name != null) {
					print("S", "s"+id, seq, "SN:Z:"+name, "SS:i:0", "SR:i:"+rank);
					++id;
				}
				name = m[1], seq.length = 0;
			} else throw Error("Wrong FASTA format!");
		}
	}
	if (name != null)
		print("S", "s"+id, seq, "SN:Z:"+name, "SS:i:0", "SR:i:"+rank);
	file.close();
	buf.destroy();
}

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
				print(">" + prefix + sep + m[1]);
			} else throw Error("Wrong FASTA format!");
		}
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
		print("  fa2gfa       convert FASTA to rGFA");
		print("  renamefa     add a prefix to sequence names in FASTA");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'fa2gfa') mg_cmd_fa2gfa(args);
	else if (cmd == 'renamefa') mg_cmd_renamefa(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
