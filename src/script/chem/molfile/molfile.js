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

	if (!this.v3000 && molecule.rgroups.size > 0) {
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
		this.writeCTab3000(this.molecule.getScaffold());
		this.writeRGroups3000(this.molecule);
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
	this.writeCollectionBlock3000(molecule);
	this.writeSGroupBlock3000(molecule);
	this.writeCR('M  V30 END CTAB');
};

Molfile.prototype.writeSGroupBlock3000 = function (molecule) {
	if (molecule.sgroups.size === 0) {
		return;
	}
	this.writeCR('M  V30 BEGIN SGROUP');
	var sGroupMap = {};
	var sgmapback = {};
	var sGroupOrder = molecule.sGroupForest.getSGroupsBFS();
	sGroupOrder.forEach((id, index) => {
		sgmapback[index + 1] = id;
		sGroupMap[id] = index + 1;
	});
	for (var q = 1; q <= sGroupOrder.size; ++q) { // each group on its own
		var id = sgmapback[q];
		var sgroup = molecule.sgroups.get(id);
		let sgDetails = `${this.getMol3000Prefix()}${this.fromKetcherIndex(id)} ${sgroup.type} 0`;
		if (sgroup.atoms.length) {
			// atoms that define the Sgroup
			// ATOMS=(natoms aidx1 aidx2 ..)
			sgDetails += ` ATOMS=(${sgroup.atoms.length} ${sgroup.atoms.map(atom => this.fromKetcherIndex(atom)).join(' ')})`;
		}
		if (sgroup.bonds && sgroup.bonds.length) {
			// crossing bonds
			// XBONDS=(nxbonds bidx1 bidx2 ..)
			sgDetails += ` XBONDS=(${sgroup.bonds.length} ${sgroup.bonds.map(bond => this.fromKetcherIndex(bond)).join(' ')})`;
		}
		if (sgroup.patoms) {
			// paradigmatic repeating unit atoms
			// PATOMS=(npatoms aidx1 aidx2 ..)
			sgDetails += ` PATOMS=(${sgroup.patoms.length} ${sgroup.patoms.map(patom => this.fromKetcherIndex(patom)).join(' ')})`;
		}
		// TODO: Add support for following:
		// 1. CBONDS
		// 2. SUBTYPE
		// 3. MULT
		// 4. COMPNO
		// 5. XBHEAD
		// 6. XBCORR
		// 7. ESTATE
		// 8. CSTATE
		// 9. CLASS
		// 10. SAP
		// 11. BRKTYP
		if (sgroup.type === 'DAT') {
			this.writeSGroupDataAttributes3000(sgroup);
		}
		sgDetails += ` CONNECT=${sgroup.data.connectivity.toUpperCase()}`;
		var parentId = this.fromKetcherIndex(this.molecule.sGroupForest.parent.get(id));
		if (parentId > 0) {
			// parent Sgroup index.
			sgDetails += ` PARENT=${parentId}`;
		}
		const brackedCoords = sgroup.bracketBox;
		if (brackedCoords) {
			// display coordinates in each bracket.
			const bracket1X = this.getAtomCoordinate3000(brackedCoords.p0.x, 4);
			const bracket1Y = this.getAtomCoordinate3000(brackedCoords.p0.y, 4);
			const bracket1Z = 0;
			const bracket2X = this.getAtomCoordinate3000(brackedCoords.p1.x, 4);
			const bracket2Y = this.getAtomCoordinate3000(brackedCoords.p1.y, 4);
			const bracket2Z = 0;
			sgDetails += ` BRKXYZ=(9 ${bracket1X} ${bracket1Y} ${bracket1Z} ${bracket2X} ${bracket2Y} ${bracket2Z} 0 0 0)`;
		}
		// display label for this Sgroup.
		sgDetails += ` LABEL=${sgroup.data.subscript}`;
		this.splitAndWriteMolV3000Lines(sgDetails);
	}
	this.writeCR('M  V30 END SGROUP');
};

/**
 * 1. Splits a single string of MOL V3000 line, into multiple lines
 *    each having length no greater than ~ 80 characters.
 * 2. Writes these strings (lines) to the MOL V3000 representation.
 *
 * From page 80 of the spec,
 *
 * To allow continuation when the 80-character line is too short, use a dash (-) as the last character.
 * When read, the line is concatenated with the next line by removing the dash
 * and stripping the initial "M V30" from the following line. For example:
 * M V30 10 20 30 “abc-
 * M V30 def”
 * is read as:
 * M V30 10 20 30 "abc def"
 *
 * @param text
 */
Molfile.prototype.splitAndWriteMolV3000Lines = function (text) {
	const textLines = text.match(/.{1,70}/g);
	this.writeCR(textLines.join('-\nM  V30 '));
};

Molfile.prototype.writeSGroupDataAttributes3000 = function (sgroup) {
	let sgDetails = '';
	if (sgroup.data.fieldName.length > 0) {
		// the name of data field for Data Sgroup.
		sgDetails += ` FIELDNAME=${sgroup.data.fieldName}`;
	}
	// TODO: confirm whether FIELDINFO is just units.
	if (sgroup.data.units && sgroup.data.units.length > 0) {
		sgDetails += ` FIELDINFO=${sgroup.data.units}`;
	}
	if (sgroup.pp) {
		// NOTE: please refer to https://drive.google.com/file/d/0Bx3dsPc7eyZKdXlHTktQTFJUMHc/view
		//
		// *Data Sgroup Display Information [Sgroup]*
		// M SDD sss xxxxx.xxxxyyyyy.yyyy eeefgh i jjjkkk ll m noo
		//
		// sss: Index of data Sgroup
		// x,y: Coordinates (2F10.4)
		// eee: (Reserved for future use)
		// f: Data display, A = attached, D = detached
		// g: Absolute, relative placement, A = absolute, R = relative
		// h: Display units, blank = no units displayed, U = display units
		// i: (Reserved for future use)
		// jjj: Number of characters to display (1...999 or ALL)
		// kkk: Number of lines to display (unused, always 1)
		// ll: (Reserved for future use)
		// m: Tag character for tagged detached display (if non-blank)
		// n: Data display DASP position (1...9). (MACCS-II only)
		// oo: (Reserved for future use)
		const x = this.getAtomCoordinate3000(sgroup.pp.x, 4).toString().padStart(10, ' ');
		const y = this.getAtomCoordinate3000(-sgroup.pp.y, 4).toString().padStart(10, ' ');
		const eee = '   ';
		const f = sgroup.data.attached ? 'A' : 'D';
		const g = sgroup.data.absolute ? 'A' : 'R';
		const h = sgroup.data.showUnits ? 'U' : ' ';
		const i = ' ';
		const jjj = sgroup.data.nCharsToDisplay === -1 ? 'ALL' : sgroup.data.nCharsToDisplay.padStart(3, 0);
		const kkk = '001';
		const ll = ' ';
		const m = sgroup.data.tagChar;
		const n = sgroup.data.daspPos;
		const oo = '  ';
		sgDetails += ` FIELDDISP="${x}${y} ${eee}${f}${g}${h} ${i} ${jjj}${kkk} ${ll} ${m}  ${n}${oo}"`;
	}
	if (sgroup.data.query.length > 0) {
		// querytype is the type of query or no query if missing.
		sgDetails += ` QUERYTYPE=${sgroup.data.query}`;
	}
	if (sgroup.data.queryOp.length > 0) {
		// queryop is the query operator
		sgDetails += ` QUERYOP=${sgroup.data.queryOp}`;
	}
	if (sgroup.data.fieldValue.length > 0) {
		// fielddata is the query or field data.
		sgDetails += ` FIELDDATA=${sgroup.data.fieldValue}`;
	}
	return sgDetails;
};


/**
 * Writes R-Group block for each R-Group in the structure
 * @param molecule
 */
Molfile.prototype.writeRGroups3000 = function (molecule) { // eslint-disable-line max-statements
	/* saver */
	if (molecule.rgroups.size === 0) {
		return;
	}
	molecule.rgroups.forEach((rgroup, id) => {
		this.writeRGroup3000(molecule, rgroup, id);
	});
};

/**
 * Writes R Group Block in V3000 format
 *
 *
 * Example of an R Group Block
 *
	 M  V30 BEGIN RGROUP 1
	 M  V30 RLOGIC 0 0 ""
	 M  V30 BEGIN CTAB
	 M  V30 COUNTS 5 5 0 0 0
	 M  V30 BEGIN ATOM
	 M  V30 1 C -2.9792 2.0183 0 0
	 M  V30 2 C -4.225 1.113 0 0
	 M  V30 3 C -3.7492 -0.3516 0 0
	 M  V30 4 C -2.2092 -0.3516 0 0
	 M  V30 5 C -1.7334 1.113 0 0
	 M  V30 END ATOM
	 M  V30 BEGIN BOND
	 M  V30 1 1 1 2
	 M  V30 2 1 1 5
	 M  V30 3 1 2 3
	 M  V30 4 1 3 4
	 M  V30 5 1 4 5
	 M  V30 END BOND
	 M  V30 END CTAB
	 M  V30 END RGROUP

 * @param molecule
 * @param rgroup
 * @param id
 */
Molfile.prototype.writeRGroup3000 = function (molecule, rgroup, id) { // eslint-disable-line max-statements
	/* saver */
	this.writeCR(`M  V30 BEGIN RGROUP ${id}`);
	this.writeCR(`M  V30 RLOGIC ${rgroup.ifthen} ${rgroup.resth ? 1 : 0} ${rgroup.range}`);
	rgroup.frags.forEach((fid) => {
		this.writeCTab3000(molecule.getFragment(fid));
	});
	this.writeCR('M  V30 END RGROUP');
};

/**
 * Writes Atom Block for V3000 format
 *
 * Example Atom Block:
 *
	 M  V30 BEGIN ATOM
	 M  V30 1 C -7.3393 7.2781 0 0
	 M  V30 2 C -8.673 6.5081 0 0
	 M  V30 3 C -6.0056 6.5081 0 0
	 M  V30 END ATOM
 * @param molecule
 */
Molfile.prototype.writeAtomBlock3000 = function (molecule) {
	if (molecule.atoms.size === 0) {
		return;
	}
	this.writeCR('M  V30 BEGIN ATOM');
	molecule.atoms.forEach((atom, id) => {
		let atomDetails = this.getMol3000Prefix();
		atomDetails += `${this.fromKetcherIndex(id)} `;
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
		if (atom.isotope !== 0) {
			atomDetails += ` MASS=${atom.isotope}`;
		}
		// TODO: Add support for following:
		// 1. CFG - Stereo configuration
		// 3. STBOX - Stereo box
		// 4. ATTCHORD - Attachment order

		this.splitAndWriteMolV3000Lines(atomDetails);
	});
	this.writeCR('M  V30 END ATOM');
};

/**
 * Ketcher's struct object maintains indices as 0-based,
 * but we write them as 1-based indices to the outside world.
 * This helper simply converts a 0-based index to its corresponding 1-based index.
 *
 * @param index
 */
Molfile.prototype.fromKetcherIndex = function (index) {
	return index + 1;
}

/**
 * Writes Bond Block for V3000 format
 *
 * Example Bond Block:
 *
	 M  V30 BEGIN BOND
	 M  V30 1 1 1 2
	 M  V30 2 1 1 3
	 M  V30 END BOND
 * @param molecule
 */
Molfile.prototype.writeBondBlock3000 = function (molecule) {
	if (molecule.bonds.size === 0) {
		return;
	}
	this.writeCR('M  V30 BEGIN BOND');
	molecule.bonds.forEach((bond, id) => {
		let bondDetails = this.getMol3000Prefix();
		bondDetails += `${this.fromKetcherIndex(id)} `;
		bondDetails += bond.type;
		bondDetails += ' ';
		bondDetails += this.fromKetcherIndex(bond.begin);
		bondDetails += ' ';
		bondDetails += this.fromKetcherIndex(bond.end);

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
		this.splitAndWriteMolV3000Lines(bondDetails);
	});
	this.writeCR('M  V30 END BOND');
};

/**
 * @param {Struct} molecule
 *
 * Example collection block:
 * 	M  V30 BEGIN COLLECTION
 * 	M  V30 MDLV30/STEABS ATOMS=(2 3 4)
 * 	M  V30 MDLV30/STEREL4 ATOMS=(1 2)
 * 	M  V30 MDLV30/STERAC5 ATOMS=(2 1 6)
 * 	M  V30 END COLLECTION
 */
Molfile.prototype.writeCollectionBlock3000 = function (molecule) {
	const { abs, rac, rel } = molecule.enhancedStereo;
	this.writeCR('M  V30 BEGIN COLLECTION');

	// Absolute stereochemistry collection
	// eg: M  V30 MDLV30/STEABS ATOMS=(2 3 4)
	if (abs.length) {
		const atomIndices = abs.map(this.fromKetcherIndex).join(' ');
		let absoluteDetails = this.getMol3000Prefix();
		absoluteDetails += `MDLV30/STEABS ATOMS=(${abs.length} ${atomIndices})`;
		const absoluteDetailsLines = absoluteDetails.match(/.{1,70}/g);
		this.writeCR(absoluteDetailsLines.join(' -\nM  V30 '));
	}

	// Relative stereochemistry collection
	// eg: M  V30 MDLV30/STEREL4 ATOMS=(1 2)
	for (let [relLabel, relIndices] of rel.entries()) {
		const atomIndices = relIndices.map(this.fromKetcherIndex).join(' ');
		let absoluteDetails = this.getMol3000Prefix();
		absoluteDetails += `MDLV30/STEREL${relLabel} ATOMS=(${relIndices.length} ${atomIndices})`;
		const absoluteDetailsLines = absoluteDetails.match(/.{1,70}/g);
		this.writeCR(absoluteDetailsLines.join(' -\nM  V30 '));
	}

	// "Racemic" stereochemistry collection
	// eg: M  V30 MDLV30/STERAC5 ATOMS=(2 1 6)
	for (let [racLabel, racIndices] of rac.entries()) {
		const atomIndices = racIndices.map(this.fromKetcherIndex).join(' ');
		let absoluteDetails = this.getMol3000Prefix();
		absoluteDetails += `MDLV30/STERAC${racLabel} ATOMS=(${racIndices.length} ${atomIndices})`;
		const absoluteDetailsLines = absoluteDetails.match(/.{1,70}/g);
		this.writeCR(absoluteDetailsLines.join(' -\nM  V30 '));
	}

	this.writeCR('M  V30 END COLLECTION');
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
		this.writePaddedNumber(this.fromKetcherIndex(atomProps_list[0].id), 3);
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
				if (atom.rglabel & (1 << rgi)) rglabelList.push([id, this.fromKetcherIndex(rgi)]);
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
			this.writePaddedNumber(this.fromKetcherIndex(aid), 4);
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
