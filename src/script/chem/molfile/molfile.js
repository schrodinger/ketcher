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
import { Bond } from '../struct';
import element from './../element';

import common from './common';
import utils from './utils';

const pkg = require('../../../../package.json');

function Molfile(v3000) {
	/* reader */
	/* saver */
	this.molecule = null;
	this.molfile = null;
	this.v3000 = v3000 || false;
}

Molfile.prototype.parseCTFile = function (molfileLines, shouldReactionRelayout) {
	let ret = null;
	if (molfileLines[0].search('\\$RXN') === 0)
		ret = common.parseRxn(molfileLines, shouldReactionRelayout);
	else
		ret = common.parseMol(molfileLines);
	ret.initHalfBonds();
	ret.initNeighbors();
	ret.markFragments();
	return ret;
};

Molfile.prototype.prepareSGroups = function (skipErrors, preserveIndigoDesc) {
	var mol = this.molecule;
	var toRemove = [];
	var errors = 0;

	this.molecule.sGroupForest.getSGroupsBFS().reverse().forEach((id) => {
		var sgroup = mol.sgroups.get(id);
		var errorIgnore = false;

		try {
			common.prepareForSaving[sgroup.type](sgroup, mol);
		} catch (ex) {
			if (!skipErrors || typeof (ex.id) != 'number')
				throw new Error(`Error: ${ex.message}`);
			errorIgnore = true;
		}
		/* eslint-disable no-mixed-operators*/
		if (errorIgnore ||
			!preserveIndigoDesc && /^INDIGO_.+_DESC$/i.test(sgroup.data.fieldName)) {
		/* eslint-enable no-mixed-operators*/
			errors += errorIgnore;
			toRemove.push(sgroup.id);
		}
	}, this);
	if (errors)
		throw new Error('Warning: ' + errors + ' invalid S-groups were detected. They will be omitted.');

	for (var i = 0; i < toRemove.length; ++i)
		mol.sGroupDelete(toRemove[i]);
	return mol;
};

Molfile.prototype.getCTab = function (molecule, rgroups) {
	/* saver */
	this.molecule = molecule.clone();
	this.prepareSGroups(false, false);
	this.molfile = '';
	this.writeCTab2000(rgroups);
	return this.molfile;
};

Molfile.prototype.saveMolecule = function (molecule, skipSGroupErrors, norgroups, preserveIndigoDesc) { // eslint-disable-line max-statements
	/* saver */
	this.reaction = molecule.rxnArrows.size > 0;
	if (molecule.rxnArrows.size > 1)
		throw new Error('Reaction may not contain more than one arrow');
	this.molfile = '' + molecule.name;
	if (this.reaction) {
		if (molecule.rgroups.size > 0)
			throw new Error('Reactions with r-groups are not supported at the moment');
		var components = molecule.getComponents();

		var reactants = components.reactants;
		var products = components.products;
		return this.v3000 ? this.writeV3000Reaction(molecule, reactants, products) :
			this.writeV2000Reaction(molecule, reactants, products);
	}

	if (molecule.rgroups.size > 0) {
		if (norgroups) {
			molecule = molecule.getScaffold();
		} else {
			var scaffold = new Molfile(false).getCTab(molecule.getScaffold(), molecule.rgroups);
			this.molfile = '$MDL  REV  1\n$MOL\n$HDR\n' + molecule.name + '\n\n\n$END HDR\n';
			this.molfile += '$CTAB\n' + scaffold + '$END CTAB\n';

			molecule.rgroups.forEach((rg, rgid) => {
				this.molfile += '$RGP\n';
				this.writePaddedNumber(rgid, 3);
				this.molfile += '\n';
				rg.frags.forEach((fid) => {
					const group = new Molfile(false).getCTab(molecule.getFragment(fid));
					this.molfile += '$CTAB\n' + group + '$END CTAB\n';
				});
				this.molfile += '$END RGP\n';
			});
			this.molfile += '$END MOL\n';

			return this.molfile;
		}
	}

	this.molecule = molecule.clone();

	this.prepareSGroups(skipSGroupErrors, preserveIndigoDesc);

	if (this.v3000) {
		this.writeHeaderV3000();
		this.writeCTab3000(this.molecule);
		this.writeCR('M  END');
	} else {
		this.writeHeader();
		this.writeCTab2000();
	}

	return this.molfile;
};

/**
 * Writes the given reaction to V2000 format.
 *
 * @param molecule {Struct}
 * @param reactants {!Array<Pile<number>>}
 * @param products {!Array<Pile<number>>}
 * @returns {string}
 */
Molfile.prototype.writeV2000Reaction = function (molecule, reactants, products) {
	var all = reactants.concat(products);
	this.molfile = '$RXN\n' + molecule.name + '\n\n\n' +
		utils.paddedNum(reactants.length, 3) +
		utils.paddedNum(products.length, 3) +
		utils.paddedNum(0, 3) + '\n';
	for (var i = 0; i < all.length; ++i) {
		var saver = new Molfile(false);
		var submol = molecule.clone(all[i], null, true);
		var molfile = saver.saveMolecule(submol, false, true);
		this.molfile += '$MOL\n' + molfile;
	}
	return this.molfile;
};

/**
 * Writes the given reaction to V3000 format.
 *
 * @param molecule {Struct}
 * @param reactants {!Array<Pile<number>>}
 * @param products {!Array<Pile<number>>}
 * @returns {string}
 */
Molfile.prototype.writeV3000Reaction = function (molecule, reactants, products) {
	this.writeCR('$RXN V3000');
	this.writeCR();
	this.writeCR('      Ketcher');
	this.writeCR();
	this.writeCR(`M  V30 COUNTS ${reactants.length} ${products.length}`);
	if (reactants.length) {
		this.writeCR('M  V30 BEGIN REACTANT');
		reactants.forEach((reactant) => {
			this.writeCTab3000(molecule.clone(reactant, null, true));
		});
		this.writeCR('M  V30 END REACTANT');
	}
	if (products.length) {
		this.writeCR('M  V30 BEGIN PRODUCT');
		products.forEach((product) => {
			this.writeCTab3000(molecule.clone(product, null, true));
		});
		this.writeCR('M  V30 END PRODUCT');
	}
	this.write('M  END');
	return this.molfile;
};

Molfile.prototype.writeHeader = function () {
	/* saver */

	var date = new Date();

	this.writeCR(); // TODO: write structure name
	this.writeWhiteSpace(2);
	this.write('Ketcher');
	this.writeWhiteSpace();
	this.writeCR(((date.getMonth() + 1) + '').padStart(2) + (date.getDate() + '').padStart(2) + ((date.getFullYear() % 100) + '').padStart(2) +
	(date.getHours() + '').padStart(2) + (date.getMinutes() + '').padStart(2) + '2D 1   1.00000     0.00000     0');
	this.writeCR();
};

Molfile.prototype.writeHeaderV3000 = function () {
	/* saver */
	this.writeCR();
	this.writeCR(`  Ketcher ${pkg.version}`);
	this.writeCR();
	this.writeCR('  0  0  0     0  0            999 V3000');
};

Molfile.prototype.write = function (str) {
	/* saver */
	this.molfile += str;
};

Molfile.prototype.writeCR = function (str) {
	/* saver */
	if (arguments.length == 0)
		str = '';

	this.molfile += str + '\n';
};

Molfile.prototype.writeWhiteSpace = function (length) {
	/* saver */

	if (arguments.length == 0)
		length = 1;

	this.write(' '.repeat(Math.max(length, 0)));
};

Molfile.prototype.writePadded = function (str, width) {
	/* saver */
	this.write(str);
	this.writeWhiteSpace(width - str.length);
};

Molfile.prototype.writePaddedNumber = function (number, width) {
	/* saver */

	var str = (number - 0).toString();

	this.writeWhiteSpace(width - str.length);
	this.write(str);
};

Molfile.prototype.writePaddedFloat = function (number, width, precision) {
	/* saver */

	this.write(utils.paddedNum(number, width, precision));
};

Molfile.prototype.writeCTab2000Header = function () {
	/* saver */

	this.writePaddedNumber(this.molecule.atoms.size, 3);
	this.writePaddedNumber(this.molecule.bonds.size, 3);

	this.writePaddedNumber(0, 3);
	this.writeWhiteSpace(3);
	this.writePaddedNumber(this.molecule.isChiral ? 1 : 0, 3);
	this.writePaddedNumber(0, 3);
	this.writeWhiteSpace(12);
	this.writePaddedNumber(999, 3);
	this.writeCR(' V2000');
};


Molfile.prototype.writeCTab3000 = function (molecule) { // eslint-disable-line max-statements
	/* saver */
	this.writeCR('M  V30 BEGIN CTAB');
	const noOfAtoms = molecule.atoms.size;
	const noOfBonds = molecule.bonds.size;
	const noOfSGroups = molecule.sgroups.size;
	const isChiral = molecule.isChiral ? 1 : 0;
	// Counts Line
	this.writeCR(`M  V30 COUNTS ${noOfAtoms} ${noOfBonds} ${noOfSGroups} 0 ${isChiral}`);
	this.writeAtomBlock3000(molecule);
	this.writeBondBlock3000(molecule);
	this.writeCR('M  V30 END CTAB');
};

Molfile.prototype.writeAtomBlock3000 = function (molecule) {
	if (molecule.atoms.size === 0) {
		return;
	}
	this.writeCR('M  V30 BEGIN ATOM');
	molecule.atoms.forEach((atom, id) => {
		let atomDetails = this.getMol3000Prefix();
		atomDetails += `${id + 1} `;
		atomDetails += this.getAtomType(atom);
		atomDetails += ' ';
		// It would be uniform to keep a maximum of 4 decimal places,
		// which should be accurate enough for positioning the molecule.
		// It would also keep the representation neat enough.
		// Technically we can have any number of decimal places
		atomDetails += this.getAtomCoordinate3000(atom.pp.x, 4);
		atomDetails += ' ';
		atomDetails += this.getAtomCoordinate3000(-atom.pp.y, 4);
		atomDetails += ' ';
		atomDetails += this.getAtomCoordinate3000(atom.pp.z, 0);
		atomDetails += ' ';
		atomDetails += atom.aam;
		if (atom.charge !== 0) {
			atomDetails += ` CHG=${atom.charge}`;
		}
		if (atom.radical !== 0) {
			atomDetails += ` RAD=${atom.radical}`;
		}
		if (atom.valence !== 0) {
			atomDetails += ` VAL=${atom.valence}`;
		}
		if (atom.hCount !== 0) {
			atomDetails += ` HCOUNT=${atom.hCount}`;
		}
		if (atom.invRet !== 0) {
			atomDetails += ` INVRET=${atom.invRet}`;
		}
		if (atom.exactChangeFlag != 0) {
			atomDetails += ` EXACHG=${atom.exactChangeFlag}`;
		}
		if (atom.substitutionCount !== 0) {
			atomDetails += ` SUBST=${atom.substitutionCount}`;
		}
		if (atom.unsaturatedAtom !== 0) {
			atomDetails += ` UNSAT=${atom.unsaturatedAtom}`;
		}
		if (atom.ringBondCount !== 0) {
			atomDetails += ` RBCNT=${atom.ringBondCount}`;
		}
		if (atom.attpnt != null) {
			atomDetails += ` ATTCHPT=${atom.attpnt}`;
		}
		if (atom.rglabel != null && atom.label === 'R#') {
			atomDetails += ` RGROUPS=(1 ${atom.rglabel})`;
		}
		// TODO: Add support for following:
		// 1. CFG - Stereo configuration
		// 2. MASS - Atomic weight
		// 3. STBOX - Stereo box
		// 4. ATTCHORD - Attachment order

		const atomDetailsLines = atomDetails.match(/.{1,70}/g);
		this.writeCR(atomDetailsLines.join(' -\nM  V30 '));
	});
	this.writeCR('M  V30 END ATOM');
};

Molfile.prototype.writeBondBlock3000 = function (molecule) {
	if (molecule.bonds.size === 0) {
		return;
	}
	this.writeCR('M  V30 BEGIN BOND');
	molecule.bonds.forEach((bond, id) => {
		let bondDetails = this.getMol3000Prefix();
		bondDetails += `${id + 1} `;
		bondDetails += bond.type;
		bondDetails += ' ';
		bondDetails += bond.begin + 1;
		bondDetails += ' ';
		bondDetails += bond.end + 1;

		if (bond.topology !== 0) {
			bondDetails += ` TOPO=${bond.topology}`;
		}
		if (bond.reactingCenterStatus !== 0) {
			bondDetails += ` RXCTR=${bond.reactingCenterStatus}`;
		}
		if (bond.stereo !== 0) {
			bondDetails += ` CFG=${this.getV3000BondConfiguration(bond.stereo)}`;
		}
		// TODO: Add support for following:
		// 1. STBOX - Stereo box
		const bondDetailsLines = bondDetails.match(/.{1,70}/g);
		this.writeCR(bondDetailsLines.join(' -\nM  V30 '));
	});
	this.writeCR('M  V30 END BOND');
};

Molfile.prototype.getMol3000Prefix = function () {
	return 'M  V30 ';
};

Molfile.prototype.getV3000BondConfiguration = function (stereo) {
	switch(stereo) {
		case Bond.PATTERN.STEREO.UP:
			return 1;
		case Bond.PATTERN.STEREO.EITHER:
			return 2;
		case Bond.PATTERN.STEREO.DOWN:
			return 3;
		default:
			return 0;
	}
};

Molfile.prototype.getAtomCoordinate3000 = function (number, precision) {
	number = parseFloat(number);
	return number.toFixed(precision || 0).replace(',', '.');
};

/**
 * Returns the atom type
 *
 * @param atom
 * @returns {string}
 */
Molfile.prototype.getAtomType = function (atom) {
	var label = atom.label;
	if (atom.atomList != null) {
		return 'L';
	}
	if (atom['pseudo']) {
		if (atom['pseudo'].length > 3) {
			return 'A';
		}
	}
	if (!element.map[label] && ['A', 'Q', 'X', '*', 'R#'].indexOf(label) == -1) { // search in generics?
		return 'C';
	}
	return label;
};

Molfile.prototype.writeCTab2000 = function (rgroups) { // eslint-disable-line max-statements
	/* saver */
	this.writeCTab2000Header();

	this.mapping = {};
	var i = 1;

	/* eslint-disable camelcase*/
	var atomList_list = [];
	var atomProps_list = [];
	/* eslint-enable camel-case*/
	this.molecule.atoms.forEach((atom, id) => {
		this.writePaddedFloat(atom.pp.x, 10, 4);
		this.writePaddedFloat(-atom.pp.y, 10, 4);
		this.writePaddedFloat(atom.pp.z, 10, 4);
		this.writeWhiteSpace();

		var label = atom.label;
		if (atom.atomList != null) {
			label = 'L';
			atomList_list.push(id);
		} else if (atom['pseudo']) {
			if (atom['pseudo'].length > 3) {
				label = 'A';
				atomProps_list.push({ id, value: '\'' + atom['pseudo'] + '\'' });
			}
		} else if (atom['alias']) {
			atomProps_list.push({ id, value: atom['alias'] });
		} else if (!element.map[label] && ['A', 'Q', 'X', '*', 'R#'].indexOf(label) == -1) { // search in generics?
			label = 'C';
			atomProps_list.push({ id, value: atom.label });
		}
		this.writePadded(label, 3);
		this.writePaddedNumber(0, 2);
		this.writePaddedNumber(0, 3);
		this.writePaddedNumber(0, 3);

		if (typeof atom.hCount === 'undefined')
			atom.hCount = 0;
		this.writePaddedNumber(atom.hCount, 3);

		if (typeof atom.stereoCare === 'undefined')
			atom.stereoCare = 0;
		this.writePaddedNumber(atom.stereoCare, 3);

		this.writePaddedNumber(atom.explicitValence < 0 ? 0 : (atom.explicitValence == 0 ? 15 : atom.explicitValence), 3); // eslint-disable-line no-nested-ternary

		this.writePaddedNumber(0, 3);
		this.writePaddedNumber(0, 3);
		this.writePaddedNumber(0, 3);

		if (typeof atom.aam === 'undefined')
			atom.aam = 0;
		this.writePaddedNumber(atom.aam, 3);

		if (typeof atom.invRet === 'undefined')
			atom.invRet = 0;
		this.writePaddedNumber(atom.invRet, 3);

		if (typeof atom.exactChangeFlag === 'undefined')
			atom.exactChangeFlag = 0;
		this.writePaddedNumber(atom.exactChangeFlag, 3);

		this.writeCR();

		this.mapping[id] = i;
		i++;
	}, this);

	this.bondMapping = {};
	i = 1;
	this.molecule.bonds.forEach((bond, id) => {
		this.bondMapping[id] = i++;
		this.writePaddedNumber(this.mapping[bond.begin], 3);
		this.writePaddedNumber(this.mapping[bond.end], 3);
		this.writePaddedNumber(bond.type, 3);

		if (typeof bond.stereo === 'undefined')
			bond.stereo = 0;
		this.writePaddedNumber(bond.stereo, 3);

		this.writePadded(bond.xxx, 3);

		if (typeof bond.topology === 'undefined')
			bond.topology = 0;
		this.writePaddedNumber(bond.topology, 3);

		if (typeof bond.reactingCenterStatus === 'undefined')
			bond.reactingCenterStatus = 0;
		this.writePaddedNumber(bond.reactingCenterStatus, 3);

		this.writeCR();
	});

	while (atomProps_list.length > 0) {
		this.write('A  ');
		this.writePaddedNumber(atomProps_list[0].id + 1, 3);
		this.writeCR();
		this.writeCR(atomProps_list[0].value);
		atomProps_list.splice(0, 1);
	}

	var chargeList = [];
	var isotopeList = [];
	var radicalList = [];
	var rglabelList = [];
	var rglogicList = [];
	var aplabelList = [];
	var rbcountList = [];
	var unsaturatedList = [];
	var substcountList = [];

	this.molecule.atoms.forEach((atom, id) => {
		if (atom.charge != 0)
			chargeList.push([id, atom.charge]);
		if (atom.isotope != 0)
			isotopeList.push([id, atom.isotope]);
		if (atom.radical != 0)
			radicalList.push([id, atom.radical]);
		if (atom.rglabel != null && atom.label == 'R#') { // TODO need to force rglabel=null when label is not 'R#'
			for (var rgi = 0; rgi < 32; rgi++)
				if (atom.rglabel & (1 << rgi)) rglabelList.push([id, rgi + 1]);
		}
		if (atom.attpnt != null)
			aplabelList.push([id, atom.attpnt]);
		if (atom.ringBondCount != 0)
			rbcountList.push([id, atom.ringBondCount]);
		if (atom.substitutionCount != 0)
			substcountList.push([id, atom.substitutionCount]);
		if (atom.unsaturatedAtom != 0)
			unsaturatedList.push([id, atom.unsaturatedAtom]);
	});

	if (rgroups) {
		rgroups.forEach((rg, rgid) => {
			if (rg.resth || rg.ifthen > 0 || rg.range.length > 0) {
				var line = '  1 ' +
					utils.paddedNum(rgid, 3) + ' ' +
					utils.paddedNum(rg.ifthen, 3) + ' ' +
					utils.paddedNum(rg.resth ? 1 : 0, 3) + '   ' + rg.range;
				rglogicList.push(line);
			}
		});
	}

	function writeAtomPropList(propId, values) {
		while (values.length > 0) {
			var part = [];

			while (values.length > 0 && part.length < 8) {
				part.push(values[0]);
				values.splice(0, 1);
			}

			this.write(propId);
			this.writePaddedNumber(part.length, 3);

			part.forEach((value) => {
				this.writeWhiteSpace();
				this.writePaddedNumber(this.mapping[value[0]], 3);
				this.writeWhiteSpace();
				this.writePaddedNumber(value[1], 3);
			});

			this.writeCR();
		}
	}

	writeAtomPropList.call(this, 'M  CHG', chargeList);
	writeAtomPropList.call(this, 'M  ISO', isotopeList);
	writeAtomPropList.call(this, 'M  RAD', radicalList);
	writeAtomPropList.call(this, 'M  RGP', rglabelList);
	for (var j = 0; j < rglogicList.length; ++j)
		this.write('M  LOG' + rglogicList[j] + '\n');

	writeAtomPropList.call(this, 'M  APO', aplabelList);
	writeAtomPropList.call(this, 'M  RBC', rbcountList);
	writeAtomPropList.call(this, 'M  SUB', substcountList);
	writeAtomPropList.call(this, 'M  UNS', unsaturatedList);

	if (atomList_list.length > 0) {
		for (j = 0; j < atomList_list.length; ++j) {
			var aid = atomList_list[j];
			var atomList = this.molecule.atoms.get(aid).atomList;
			this.write('M  ALS');
			this.writePaddedNumber(aid + 1, 4);
			this.writePaddedNumber(atomList.ids.length, 3);
			this.writeWhiteSpace();
			this.write(atomList.notList ? 'T' : 'F');

			var labelList = atomList.labelList();
			for (var k = 0; k < labelList.length; ++k) {
				this.writeWhiteSpace();
				this.writePadded(labelList[k], 3);
			}
			this.writeCR();
		}
	}

	var sgmap = {};
	var cnt = 1;
	var sgmapback = {};
	var sgorder = this.molecule.sGroupForest.getSGroupsBFS();
	sgorder.forEach((id) => {
		sgmapback[cnt] = id;
		sgmap[id] = cnt++;
	});
	for (var q = 1; q < cnt; ++q) { // each group on its own
		var id = sgmapback[q];
		var sgroup = this.molecule.sgroups.get(id);
		this.write('M  STY');
		this.writePaddedNumber(1, 3);
		this.writeWhiteSpace(1);
		this.writePaddedNumber(q, 3);
		this.writeWhiteSpace(1);
		this.writePadded(sgroup.type, 3);
		this.writeCR();

		// TODO: write subtype, M SST

		this.write('M  SLB');
		this.writePaddedNumber(1, 3);
		this.writeWhiteSpace(1);
		this.writePaddedNumber(q, 3);
		this.writeWhiteSpace(1);
		this.writePaddedNumber(q, 3);
		this.writeCR();

		var parentid = this.molecule.sGroupForest.parent.get(id);
		if (parentid >= 0) {
			this.write('M  SPL');
			this.writePaddedNumber(1, 3);
			this.writeWhiteSpace(1);
			this.writePaddedNumber(q, 3);
			this.writeWhiteSpace(1);
			this.writePaddedNumber(sgmap[parentid], 3);
			this.writeCR();
		}

		// connectivity
		if (sgroup.type == 'SRU' && sgroup.data.connectivity) {
			var connectivity = '';
			connectivity += ' ';
			connectivity += q.toString().padStart(3);
			connectivity += ' ';
			connectivity += (sgroup.data.connectivity || '').padEnd(3);
			this.write('M  SCN');
			this.writePaddedNumber(1, 3);
			this.write(connectivity.toUpperCase());
			this.writeCR();
		}

		if (sgroup.type == 'SRU') {
			this.write('M  SMT ');
			this.writePaddedNumber(q, 3);
			this.writeWhiteSpace();
			this.write(sgroup.data.subscript || 'n');
			this.writeCR();
		}

		this.writeCR(common.saveToMolfile[sgroup.type](sgroup, this.molecule, sgmap, this.mapping, this.bondMapping));
	}

	// TODO: write M  APO
	// TODO: write M  AAL
	// TODO: write M  RGP
	// TODO: write M  LOG

	this.writeCR('M  END');
};

export default Molfile;
