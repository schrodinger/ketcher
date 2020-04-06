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

import Pool from '../../util/pool';

function EnhancedStereo(molecule) {
	this.abs = []; // absolute stereo info
	this.rac = new Pool(); // n -> list of atom ids // racemic stereo info
	this.rel = new Pool(); // n -> list of atom ids // relative stereo info
	this.molecule = molecule;
}

/**
 * @return {EnhancedStereo}
 */
EnhancedStereo.prototype.clone = function() {
	const cp = new EnhancedStereo(this.molecule);
	cp.abs = this.abs.slice();
	cp.rac = new Pool(this.rac);
	cp.rel = new Pool(this.rel);
	return cp;
};

export default EnhancedStereo;
