/****************************************************************************
 * Copyright 2018 EPAM Systems
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ***************************************************************************/

import 'babel-polyfill';
import 'whatwg-fetch';
import queryString from 'query-string';

import api from './api';
import molfile from './chem/molfile';
import smiles from './chem/smiles';
import toolbarTypeMap from './ui/state/toolbar/toolbartypes';
import * as structformat from './ui/data/convert/structformat';

import ui from './ui';
import Render from './render';

function getSmiles() {
	return smiles.stringify(ketcher.editor.struct(),
		{ ignoreErrors: true });
}

function saveSmiles() {
	const struct = ketcher.editor.struct();
	return structformat.toString(struct, 'smiles-ext', ketcher.server)
		.catch(() => smiles.stringify(struct));
}

/**
 * @param {boolean} v3000 - whether to retrieve the molecule in v3000 format. Default is v2000
 * @returns {string}
 */
function getMolfile(v3000) {
	return molfile.stringify(ketcher.editor.struct(),
		{ ignoreErrors: true, v3000 });
}

/**
 * Exports representation of the current sketch in Ketcher, in the specified format.
 * @param format
 * @returns {Promise<string>}
 */
function getRepresentationInFormat(format) {
	const input_mol = molfile.stringify(ketcher.editor.struct(), {
		ignoreErrors: true
	});
	// TODO: Remove Plexus usage for converting between formats
	return new Promise((resolve) => {
		$.ajax({
			type: "POST",
			url: '/plexus/rest-v0/util/calculate/stringMolExport',
			data: { structure: input_mol, parameters: format }
		}).done(function updateInputCompound(representation) {
			resolve(representation);
		}).fail(function onError() {
			throw new Error('Could not perform the conversion');
		});
	});
}

function getSvg() {
	const mol = Module.get_mol(molfile.stringify(ketcher.editor.struct(), {
		ignoreErrors: true
	}));
	return mol.get_svg();
}

function setMolecule(molString) {
	if (!(typeof molString === 'string')) {
		return;
	}
	$.ajax({
		type: "POST",
		url: '/plexus/rest-v0/util/calculate/stringMolExport',
		data: { structure: molString, parameters: 'mol:V3' }
	}).done(function updateInputCompound(molfile) {
		ketcher.ui.load(molfile, {
			rescale: true
		});
	}).fail(function onError() {
		ketcher.ui.load(molString, {
			rescale: true
		});
	});
}

function addFragment(molString) {
	if (!(typeof molString === 'string'))
		return;
	ketcher.ui.load(molString, {
		rescale: true,
		fragment: true
	});
}

function showMolfile(clientArea, molString, options) {
	const render = new Render(clientArea, Object.assign({
		scale: options.bondLength || 75
	}, options));
	if (molString) {
		const mol = molfile.parse(molString);
		render.setMolecule(mol);
	}
	render.update();
	// not sure we need to expose guts
	return render;
}

/**
 * Sets the toolbar type for Ketcher
 * @param toolbarType {toolbarTypeMap}
 */
function setToolbarType(toolbarType) {
	ketcher.ui.changeToolbarType(toolbarType);
}

// TODO: replace window.onload with something like <https://github.com/ded/domready>
// to start early
window.onload = function () {
	const params = queryString.parse(document.location.search);
	if (params.api_path)
		ketcher.apiPath = params.api_path;
	ketcher.server = api(ketcher.apiPath, {
		'smart-layout': true,
		'ignore-stereochemistry-errors': true,
		'mass-skip-error-on-pseudoatoms': false,
		'gross-formula-add-rsites': true
	});
	ketcher.ui = ui(Object.assign({}, params, buildInfo), ketcher.server);
	ketcher.editor = global._ui_editor;
	ketcher.server.then(() => {
		if (params.mol)
			ketcher.ui.load(params.mol);
	}, () => {
		document.title += ' (standalone)';
	});
};

const buildInfo = {
	version: '__VERSION__',
	apiPath: '__API_PATH__',
	buildDate: '__BUILD_DATE__',
	buildNumber: '__BUILD_NUMBER__' || null
};

const ketcher = module.exports = Object.assign({ // eslint-disable-line no-multi-assign
	getSmiles,
	getSvg,
	saveSmiles,
	getMolfile,
	getRepresentationInFormat,
	setMolecule,
	setToolbarType,
	addFragment,
	showMolfile
}, buildInfo);
