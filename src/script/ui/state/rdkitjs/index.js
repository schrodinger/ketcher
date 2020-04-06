import molfile from '../../../chem/molfile';
import { load } from '../shared';
import Action from '../../../editor/shared/action';
import { fromSgroupAddition, fromSgroupDeletion } from '../../../editor/actions/sgroup';
import Vec2 from '../../../util/vec2';

/**
 * Performs transformation on the Molecule present in Ketcher, using utilities
 * from the RDKit JS library.
 *
 * @param method: the RdkitJS transformation method to call
 * @returns {Function}
 */
export function rdkitJsTransform(method) {
	return (dispatch, getState) => {
		const state = getState();
		const struct = state.editor.struct().clone(null, null, false, new Map());
		const hasExplicitHydrogen = struct.hasExplicitHydrogen();
		const mol = Module.get_mol(molfile.stringify(struct, // eslint-disable-line no-undef
			{ ignoreErrors: true }));
		if (!mol.is_valid()) {
			alert(`The specified input molecule is invalid.`);
			return;
		}
		dispatch(load(doTransformation(method, mol, hasExplicitHydrogen), {
			rescale: method === 'layout',
			reactionRelayout: method === 'clean'
		}));
	};
}


export function rdkitJsCIP() {
	return (dispatch, getState) => {
		const state = getState();
		const struct = state.editor.struct().clone(null, null, false, new Map());
		const mol = Module.get_mol(molfile.stringify(struct, // eslint-disable-line no-undef
			{ ignoreErrors: true }));
		if (!mol.is_valid()) {
			alert(`The specified input molecule is invalid.`);
			return;
		}
		const res = JSON.parse(mol.get_stereo_tags());
		dispatch(calculateCip(res));
	};
}

export function calculateCip(result) {
	return (dispatch, getState) => {
		const state = getState();
		const editor = state.editor;
		const restruct = editor.render.ctab;
		const atomsIterator = restruct.molecule.atoms.keys();
		const keys = [...atomsIterator];
		const attributes = {
			absolute: false,
			attached: false,
			context: 'Atom',
			fieldName: 'CIP_DESC',
			fieldValue: '(R)',
			init: true
		};
		const action = new Action();
		deleteAllSGroupsWithName(restruct, action, attributes.fieldName);
		result.CIP_atoms.forEach((a) => {
			const [index, fieldValue] = a;
			const atomKey = keys[index];
			const atom = restruct.molecule.atoms.get(atomKey);
			attributes.fieldValue = fieldValue;
			action.mergeWith(fromSgroupAddition(restruct, 'DAT', [atomIdx], attributes, undefined, atom.pp));
		});
		result.CIP_bonds.forEach((b) => {
			const [atom1Index, atom2Index] = b;
			const atom1Idx = keys[atom1Index];
			const atom2Idx = keys[atom2Index];
			const atom1 = restruct.molecule.atoms.get(atom1Idx);
			const atom2 = restruct.molecule.atoms.get(atom2Idx);
			attributes.fieldValue = b[2];
			const pp = new Vec2((atom1.pp.x + atom2.pp.x) * 0.5,
				(atom1.pp.y + atom2.pp.y) * 0.5);
			action.mergeWith(fromSgroupAddition(restruct, 'DAT', [atom1Idx, atom2Idx], attributes, undefined, pp));
		});
		editor.update(action);
	};
}

function deleteAllSGroupsWithName(restruct, action, fieldName) {
	restruct.molecule.sgroups.forEach((sg, id) => {
		if (sg.data.fieldName === fieldName)
			action.mergeWith(fromSgroupDeletion(restruct, id));
	});
}

/**
 * @param method {string}
 * @param mol {*} TODO (pradeep): Figure out a typedef for mol
 * @param hasExplicitHydrogen {boolean}
 * @return {*}
 */
function doTransformation(method, mol, hasExplicitHydrogen) {
	switch (method) {
		case 'aromatize':
			return mol.get_aromatic_form();
		case 'dearomatize':
			return mol.get_kekule_form();
		case 'clean':
			return mol.get_new_coords();
		case 'toggleExplicitHydrogen':
			if (hasExplicitHydrogen) {
				return mol.remove_hs();
			} else {
				return mol.add_hs();
			}
		default:
			return '';
	}
}
