import { rdkitJsTransform, rdkitJsCIP } from '../state/rdkitjs';

export default {
	arom: {
		title: 'Aromatize',
		action: {
			thunk: rdkitJsTransform('aromatize')
		},
		className: 'rdkitjs',
		disabled: () => true
	},
	dearom: {
		title: 'Dearomatize',
		action: {
			thunk: rdkitJsTransform('dearomatize')
		},
		className: 'rdkitjs',
		disabled: () => true
	},
	clean: {
		shortcut: 'Mod+Shift+l',
		title: 'Clean Up',
		action: {
			thunk: rdkitJsTransform('clean')
		},
		className: 'rdkitjs',
		disabled: () => true
	},
	cip: {
		shortcut: 'Mod+p',
		title: 'Calculate CIP',
		action: {
			thunk: rdkitJsCIP()
		},
		className: 'rdkitjs',
		disabled: () => true
	},
	explicitHydrogen: {
		title: 'Add/Remove Explicit Hydrogen',
		action: {
			thunk: rdkitJsTransform('toggleExplicitHydrogen')
		},
		className: 'rdkitjs',
		disabled: () => true
	}
};
