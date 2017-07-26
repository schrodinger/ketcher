#!/usr/bin/env nodejs
/* eslint-env node */

var path = require('path');
var cp = require('child_process');

var minimist = require('minimist');

var opts = process.argv.slice(3);
var script = require.resolve(path.join(__dirname, process.argv[2]));
var binPath = cp.execSync('npm bin').toString().trim();
var tapScript = path.join(binPath, 'tap');

var filterSpec = {
	string: ['reporter'],
	default: { reporter: 'classic' },
	alias: { reporter: 'R' },
	unknown: o => filterSpec.string.includes(o)
};
var filterOpts = minimist(process.argv.slice(3), filterSpec);

cp.spawn(`node ${script} ${opts.join(' ')} |` +
         `${tapScript} --reporter=${filterOpts.reporter} -`,
         { shell: true, stdio: 'inherit' });