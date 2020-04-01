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

import { serverTransform } from '../state/server';

export default {
	layout: {
		shortcut: 'Mod+l',
		title: 'Layout',
		action: {
			thunk: serverTransform('layout')
		},
		disabled: (editor, server, options) => !options.app.server,
	},
	check: {
		title: 'Check Structure',
		action: { dialog: 'check' },
		disabled: (editor, server, options) => !options.app.server
	},
	analyse: {
		title: 'Calculated Values',
		action: { dialog: 'analyse' },
		disabled: (editor, server, options) => !options.app.server
	},
	recognize: {
		title: 'Recognize Molecule',
		action: { dialog: 'recognize' },
		disabled: (editor, server, options) => !options.app.server
	},
	miew: {
		title: '3D Viewer',
		action: { dialog: 'miew' },
		disabled: () => !window.Miew
	}
};
