<!DOCTYPE html>
<!-- saved from url=(0081)http://localhost:8888/nbconvert/html/mdm_model_formula_check.ipynb?download=false -->
<html><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>mdm_model_formula_check</title><script src="./mdm_model_formula_check_files/require.min.js.download"></script>
<script src="./mdm_model_formula_check_files/jquery.min.js.download"></script>

<style type="text/css">
    /*!
*
* Twitter Bootstrap
*
*/
/*!
 * Bootstrap v3.3.7 (http://getbootstrap.com)
 * Copyright 2011-2016 Twitter, Inc.
 * Licensed under MIT (https://github.com/twbs/bootstrap/blob/master/LICENSE)
 */
/*! normalize.css v3.0.3 | MIT License | github.com/necolas/normalize.css */
html {
  font-family: sans-serif;
  -ms-text-size-adjust: 100%;
  -webkit-text-size-adjust: 100%;
}
body {
  margin: 0;
}
article,
aside,
details,
figcaption,
figure,
footer,
header,
hgroup,
main,
menu,
nav,
section,
summary {
  display: block;
}
audio,
canvas,
progress,
video {
  display: inline-block;
  vertical-align: baseline;
}
audio:not([controls]) {
  display: none;
  height: 0;
}
[hidden],
template {
  display: none;
}
a {
  background-color: transparent;
}
a:active,
a:hover {
  outline: 0;
}
abbr[title] {
  border-bottom: 1px dotted;
}
b,
strong {
  font-weight: bold;
}
dfn {
  font-style: italic;
}
h1 {
  font-size: 2em;
  margin: 0.67em 0;
}
mark {
  background: #ff0;
  color: #000;
}
small {
  font-size: 80%;
}
sub,
sup {
  font-size: 75%;
  line-height: 0;
  position: relative;
  vertical-align: baseline;
}
sup {
  top: -0.5em;
}
sub {
  bottom: -0.25em;
}
img {
  border: 0;
}
svg:not(:root) {
  overflow: hidden;
}
figure {
  margin: 1em 40px;
}
hr {
  box-sizing: content-box;
  height: 0;
}
pre {
  overflow: auto;
}
code,
kbd,
pre,
samp {
  font-family: monospace, monospace;
  font-size: 1em;
}
button,
input,
optgroup,
select,
textarea {
  color: inherit;
  font: inherit;
  margin: 0;
}
button {
  overflow: visible;
}
button,
select {
  text-transform: none;
}
button,
html input[type="button"],
input[type="reset"],
input[type="submit"] {
  -webkit-appearance: button;
  cursor: pointer;
}
button[disabled],
html input[disabled] {
  cursor: default;
}
button::-moz-focus-inner,
input::-moz-focus-inner {
  border: 0;
  padding: 0;
}
input {
  line-height: normal;
}
input[type="checkbox"],
input[type="radio"] {
  box-sizing: border-box;
  padding: 0;
}
input[type="number"]::-webkit-inner-spin-button,
input[type="number"]::-webkit-outer-spin-button {
  height: auto;
}
input[type="search"] {
  -webkit-appearance: textfield;
  box-sizing: content-box;
}
input[type="search"]::-webkit-search-cancel-button,
input[type="search"]::-webkit-search-decoration {
  -webkit-appearance: none;
}
fieldset {
  border: 1px solid #c0c0c0;
  margin: 0 2px;
  padding: 0.35em 0.625em 0.75em;
}
legend {
  border: 0;
  padding: 0;
}
textarea {
  overflow: auto;
}
optgroup {
  font-weight: bold;
}
table {
  border-collapse: collapse;
  border-spacing: 0;
}
td,
th {
  padding: 0;
}
/*! Source: https://github.com/h5bp/html5-boilerplate/blob/master/src/css/main.css */
@media print {
  *,
  *:before,
  *:after {
    background: transparent !important;
    color: #000 !important;
    box-shadow: none !important;
    text-shadow: none !important;
  }
  a,
  a:visited {
    text-decoration: underline;
  }
  a[href]:after {
    content: " (" attr(href) ")";
  }
  abbr[title]:after {
    content: " (" attr(title) ")";
  }
  a[href^="#"]:after,
  a[href^="javascript:"]:after {
    content: "";
  }
  pre,
  blockquote {
    border: 1px solid #999;
    page-break-inside: avoid;
  }
  thead {
    display: table-header-group;
  }
  tr,
  img {
    page-break-inside: avoid;
  }
  img {
    max-width: 100% !important;
  }
  p,
  h2,
  h3 {
    orphans: 3;
    widows: 3;
  }
  h2,
  h3 {
    page-break-after: avoid;
  }
  .navbar {
    display: none;
  }
  .btn > .caret,
  .dropup > .btn > .caret {
    border-top-color: #000 !important;
  }
  .label {
    border: 1px solid #000;
  }
  .table {
    border-collapse: collapse !important;
  }
  .table td,
  .table th {
    background-color: #fff !important;
  }
  .table-bordered th,
  .table-bordered td {
    border: 1px solid #ddd !important;
  }
}
@font-face {
  font-family: 'Glyphicons Halflings';
  src: url('../components/bootstrap/fonts/glyphicons-halflings-regular.eot');
  src: url('../components/bootstrap/fonts/glyphicons-halflings-regular.eot?#iefix') format('embedded-opentype'), url('../components/bootstrap/fonts/glyphicons-halflings-regular.woff2') format('woff2'), url('../components/bootstrap/fonts/glyphicons-halflings-regular.woff') format('woff'), url('../components/bootstrap/fonts/glyphicons-halflings-regular.ttf') format('truetype'), url('../components/bootstrap/fonts/glyphicons-halflings-regular.svg#glyphicons_halflingsregular') format('svg');
}
.glyphicon {
  position: relative;
  top: 1px;
  display: inline-block;
  font-family: 'Glyphicons Halflings';
  font-style: normal;
  font-weight: normal;
  line-height: 1;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
.glyphicon-asterisk:before {
  content: "\002a";
}
.glyphicon-plus:before {
  content: "\002b";
}
.glyphicon-euro:before,
.glyphicon-eur:before {
  content: "\20ac";
}
.glyphicon-minus:before {
  content: "\2212";
}
.glyphicon-cloud:before {
  content: "\2601";
}
.glyphicon-envelope:before {
  content: "\2709";
}
.glyphicon-pencil:before {
  content: "\270f";
}
.glyphicon-glass:before {
  content: "\e001";
}
.glyphicon-music:before {
  content: "\e002";
}
.glyphicon-search:before {
  content: "\e003";
}
.glyphicon-heart:before {
  content: "\e005";
}
.glyphicon-star:before {
  content: "\e006";
}
.glyphicon-star-empty:before {
  content: "\e007";
}
.glyphicon-user:before {
  content: "\e008";
}
.glyphicon-film:before {
  content: "\e009";
}
.glyphicon-th-large:before {
  content: "\e010";
}
.glyphicon-th:before {
  content: "\e011";
}
.glyphicon-th-list:before {
  content: "\e012";
}
.glyphicon-ok:before {
  content: "\e013";
}
.glyphicon-remove:before {
  content: "\e014";
}
.glyphicon-zoom-in:before {
  content: "\e015";
}
.glyphicon-zoom-out:before {
  content: "\e016";
}
.glyphicon-off:before {
  content: "\e017";
}
.glyphicon-signal:before {
  content: "\e018";
}
.glyphicon-cog:before {
  content: "\e019";
}
.glyphicon-trash:before {
  content: "\e020";
}
.glyphicon-home:before {
  content: "\e021";
}
.glyphicon-file:before {
  content: "\e022";
}
.glyphicon-time:before {
  content: "\e023";
}
.glyphicon-road:before {
  content: "\e024";
}
.glyphicon-download-alt:before {
  content: "\e025";
}
.glyphicon-download:before {
  content: "\e026";
}
.glyphicon-upload:before {
  content: "\e027";
}
.glyphicon-inbox:before {
  content: "\e028";
}
.glyphicon-play-circle:before {
  content: "\e029";
}
.glyphicon-repeat:before {
  content: "\e030";
}
.glyphicon-refresh:before {
  content: "\e031";
}
.glyphicon-list-alt:before {
  content: "\e032";
}
.glyphicon-lock:before {
  content: "\e033";
}
.glyphicon-flag:before {
  content: "\e034";
}
.glyphicon-headphones:before {
  content: "\e035";
}
.glyphicon-volume-off:before {
  content: "\e036";
}
.glyphicon-volume-down:before {
  content: "\e037";
}
.glyphicon-volume-up:before {
  content: "\e038";
}
.glyphicon-qrcode:before {
  content: "\e039";
}
.glyphicon-barcode:before {
  content: "\e040";
}
.glyphicon-tag:before {
  content: "\e041";
}
.glyphicon-tags:before {
  content: "\e042";
}
.glyphicon-book:before {
  content: "\e043";
}
.glyphicon-bookmark:before {
  content: "\e044";
}
.glyphicon-print:before {
  content: "\e045";
}
.glyphicon-camera:before {
  content: "\e046";
}
.glyphicon-font:before {
  content: "\e047";
}
.glyphicon-bold:before {
  content: "\e048";
}
.glyphicon-italic:before {
  content: "\e049";
}
.glyphicon-text-height:before {
  content: "\e050";
}
.glyphicon-text-width:before {
  content: "\e051";
}
.glyphicon-align-left:before {
  content: "\e052";
}
.glyphicon-align-center:before {
  content: "\e053";
}
.glyphicon-align-right:before {
  content: "\e054";
}
.glyphicon-align-justify:before {
  content: "\e055";
}
.glyphicon-list:before {
  content: "\e056";
}
.glyphicon-indent-left:before {
  content: "\e057";
}
.glyphicon-indent-right:before {
  content: "\e058";
}
.glyphicon-facetime-video:before {
  content: "\e059";
}
.glyphicon-picture:before {
  content: "\e060";
}
.glyphicon-map-marker:before {
  content: "\e062";
}
.glyphicon-adjust:before {
  content: "\e063";
}
.glyphicon-tint:before {
  content: "\e064";
}
.glyphicon-edit:before {
  content: "\e065";
}
.glyphicon-share:before {
  content: "\e066";
}
.glyphicon-check:before {
  content: "\e067";
}
.glyphicon-move:before {
  content: "\e068";
}
.glyphicon-step-backward:before {
  content: "\e069";
}
.glyphicon-fast-backward:before {
  content: "\e070";
}
.glyphicon-backward:before {
  content: "\e071";
}
.glyphicon-play:before {
  content: "\e072";
}
.glyphicon-pause:before {
  content: "\e073";
}
.glyphicon-stop:before {
  content: "\e074";
}
.glyphicon-forward:before {
  content: "\e075";
}
.glyphicon-fast-forward:before {
  content: "\e076";
}
.glyphicon-step-forward:before {
  content: "\e077";
}
.glyphicon-eject:before {
  content: "\e078";
}
.glyphicon-chevron-left:before {
  content: "\e079";
}
.glyphicon-chevron-right:before {
  content: "\e080";
}
.glyphicon-plus-sign:before {
  content: "\e081";
}
.glyphicon-minus-sign:before {
  content: "\e082";
}
.glyphicon-remove-sign:before {
  content: "\e083";
}
.glyphicon-ok-sign:before {
  content: "\e084";
}
.glyphicon-question-sign:before {
  content: "\e085";
}
.glyphicon-info-sign:before {
  content: "\e086";
}
.glyphicon-screenshot:before {
  content: "\e087";
}
.glyphicon-remove-circle:before {
  content: "\e088";
}
.glyphicon-ok-circle:before {
  content: "\e089";
}
.glyphicon-ban-circle:before {
  content: "\e090";
}
.glyphicon-arrow-left:before {
  content: "\e091";
}
.glyphicon-arrow-right:before {
  content: "\e092";
}
.glyphicon-arrow-up:before {
  content: "\e093";
}
.glyphicon-arrow-down:before {
  content: "\e094";
}
.glyphicon-share-alt:before {
  content: "\e095";
}
.glyphicon-resize-full:before {
  content: "\e096";
}
.glyphicon-resize-small:before {
  content: "\e097";
}
.glyphicon-exclamation-sign:before {
  content: "\e101";
}
.glyphicon-gift:before {
  content: "\e102";
}
.glyphicon-leaf:before {
  content: "\e103";
}
.glyphicon-fire:before {
  content: "\e104";
}
.glyphicon-eye-open:before {
  content: "\e105";
}
.glyphicon-eye-close:before {
  content: "\e106";
}
.glyphicon-warning-sign:before {
  content: "\e107";
}
.glyphicon-plane:before {
  content: "\e108";
}
.glyphicon-calendar:before {
  content: "\e109";
}
.glyphicon-random:before {
  content: "\e110";
}
.glyphicon-comment:before {
  content: "\e111";
}
.glyphicon-magnet:before {
  content: "\e112";
}
.glyphicon-chevron-up:before {
  content: "\e113";
}
.glyphicon-chevron-down:before {
  content: "\e114";
}
.glyphicon-retweet:before {
  content: "\e115";
}
.glyphicon-shopping-cart:before {
  content: "\e116";
}
.glyphicon-folder-close:before {
  content: "\e117";
}
.glyphicon-folder-open:before {
  content: "\e118";
}
.glyphicon-resize-vertical:before {
  content: "\e119";
}
.glyphicon-resize-horizontal:before {
  content: "\e120";
}
.glyphicon-hdd:before {
  content: "\e121";
}
.glyphicon-bullhorn:before {
  content: "\e122";
}
.glyphicon-bell:before {
  content: "\e123";
}
.glyphicon-certificate:before {
  content: "\e124";
}
.glyphicon-thumbs-up:before {
  content: "\e125";
}
.glyphicon-thumbs-down:before {
  content: "\e126";
}
.glyphicon-hand-right:before {
  content: "\e127";
}
.glyphicon-hand-left:before {
  content: "\e128";
}
.glyphicon-hand-up:before {
  content: "\e129";
}
.glyphicon-hand-down:before {
  content: "\e130";
}
.glyphicon-circle-arrow-right:before {
  content: "\e131";
}
.glyphicon-circle-arrow-left:before {
  content: "\e132";
}
.glyphicon-circle-arrow-up:before {
  content: "\e133";
}
.glyphicon-circle-arrow-down:before {
  content: "\e134";
}
.glyphicon-globe:before {
  content: "\e135";
}
.glyphicon-wrench:before {
  content: "\e136";
}
.glyphicon-tasks:before {
  content: "\e137";
}
.glyphicon-filter:before {
  content: "\e138";
}
.glyphicon-briefcase:before {
  content: "\e139";
}
.glyphicon-fullscreen:before {
  content: "\e140";
}
.glyphicon-dashboard:before {
  content: "\e141";
}
.glyphicon-paperclip:before {
  content: "\e142";
}
.glyphicon-heart-empty:before {
  content: "\e143";
}
.glyphicon-link:before {
  content: "\e144";
}
.glyphicon-phone:before {
  content: "\e145";
}
.glyphicon-pushpin:before {
  content: "\e146";
}
.glyphicon-usd:before {
  content: "\e148";
}
.glyphicon-gbp:before {
  content: "\e149";
}
.glyphicon-sort:before {
  content: "\e150";
}
.glyphicon-sort-by-alphabet:before {
  content: "\e151";
}
.glyphicon-sort-by-alphabet-alt:before {
  content: "\e152";
}
.glyphicon-sort-by-order:before {
  content: "\e153";
}
.glyphicon-sort-by-order-alt:before {
  content: "\e154";
}
.glyphicon-sort-by-attributes:before {
  content: "\e155";
}
.glyphicon-sort-by-attributes-alt:before {
  content: "\e156";
}
.glyphicon-unchecked:before {
  content: "\e157";
}
.glyphicon-expand:before {
  content: "\e158";
}
.glyphicon-collapse-down:before {
  content: "\e159";
}
.glyphicon-collapse-up:before {
  content: "\e160";
}
.glyphicon-log-in:before {
  content: "\e161";
}
.glyphicon-flash:before {
  content: "\e162";
}
.glyphicon-log-out:before {
  content: "\e163";
}
.glyphicon-new-window:before {
  content: "\e164";
}
.glyphicon-record:before {
  content: "\e165";
}
.glyphicon-save:before {
  content: "\e166";
}
.glyphicon-open:before {
  content: "\e167";
}
.glyphicon-saved:before {
  content: "\e168";
}
.glyphicon-import:before {
  content: "\e169";
}
.glyphicon-export:before {
  content: "\e170";
}
.glyphicon-send:before {
  content: "\e171";
}
.glyphicon-floppy-disk:before {
  content: "\e172";
}
.glyphicon-floppy-saved:before {
  content: "\e173";
}
.glyphicon-floppy-remove:before {
  content: "\e174";
}
.glyphicon-floppy-save:before {
  content: "\e175";
}
.glyphicon-floppy-open:before {
  content: "\e176";
}
.glyphicon-credit-card:before {
  content: "\e177";
}
.glyphicon-transfer:before {
  content: "\e178";
}
.glyphicon-cutlery:before {
  content: "\e179";
}
.glyphicon-header:before {
  content: "\e180";
}
.glyphicon-compressed:before {
  content: "\e181";
}
.glyphicon-earphone:before {
  content: "\e182";
}
.glyphicon-phone-alt:before {
  content: "\e183";
}
.glyphicon-tower:before {
  content: "\e184";
}
.glyphicon-stats:before {
  content: "\e185";
}
.glyphicon-sd-video:before {
  content: "\e186";
}
.glyphicon-hd-video:before {
  content: "\e187";
}
.glyphicon-subtitles:before {
  content: "\e188";
}
.glyphicon-sound-stereo:before {
  content: "\e189";
}
.glyphicon-sound-dolby:before {
  content: "\e190";
}
.glyphicon-sound-5-1:before {
  content: "\e191";
}
.glyphicon-sound-6-1:before {
  content: "\e192";
}
.glyphicon-sound-7-1:before {
  content: "\e193";
}
.glyphicon-copyright-mark:before {
  content: "\e194";
}
.glyphicon-registration-mark:before {
  content: "\e195";
}
.glyphicon-cloud-download:before {
  content: "\e197";
}
.glyphicon-cloud-upload:before {
  content: "\e198";
}
.glyphicon-tree-conifer:before {
  content: "\e199";
}
.glyphicon-tree-deciduous:before {
  content: "\e200";
}
.glyphicon-cd:before {
  content: "\e201";
}
.glyphicon-save-file:before {
  content: "\e202";
}
.glyphicon-open-file:before {
  content: "\e203";
}
.glyphicon-level-up:before {
  content: "\e204";
}
.glyphicon-copy:before {
  content: "\e205";
}
.glyphicon-paste:before {
  content: "\e206";
}
.glyphicon-alert:before {
  content: "\e209";
}
.glyphicon-equalizer:before {
  content: "\e210";
}
.glyphicon-king:before {
  content: "\e211";
}
.glyphicon-queen:before {
  content: "\e212";
}
.glyphicon-pawn:before {
  content: "\e213";
}
.glyphicon-bishop:before {
  content: "\e214";
}
.glyphicon-knight:before {
  content: "\e215";
}
.glyphicon-baby-formula:before {
  content: "\e216";
}
.glyphicon-tent:before {
  content: "\26fa";
}
.glyphicon-blackboard:before {
  content: "\e218";
}
.glyphicon-bed:before {
  content: "\e219";
}
.glyphicon-apple:before {
  content: "\f8ff";
}
.glyphicon-erase:before {
  content: "\e221";
}
.glyphicon-hourglass:before {
  content: "\231b";
}
.glyphicon-lamp:before {
  content: "\e223";
}
.glyphicon-duplicate:before {
  content: "\e224";
}
.glyphicon-piggy-bank:before {
  content: "\e225";
}
.glyphicon-scissors:before {
  content: "\e226";
}
.glyphicon-bitcoin:before {
  content: "\e227";
}
.glyphicon-btc:before {
  content: "\e227";
}
.glyphicon-xbt:before {
  content: "\e227";
}
.glyphicon-yen:before {
  content: "\00a5";
}
.glyphicon-jpy:before {
  content: "\00a5";
}
.glyphicon-ruble:before {
  content: "\20bd";
}
.glyphicon-rub:before {
  content: "\20bd";
}
.glyphicon-scale:before {
  content: "\e230";
}
.glyphicon-ice-lolly:before {
  content: "\e231";
}
.glyphicon-ice-lolly-tasted:before {
  content: "\e232";
}
.glyphicon-education:before {
  content: "\e233";
}
.glyphicon-option-horizontal:before {
  content: "\e234";
}
.glyphicon-option-vertical:before {
  content: "\e235";
}
.glyphicon-menu-hamburger:before {
  content: "\e236";
}
.glyphicon-modal-window:before {
  content: "\e237";
}
.glyphicon-oil:before {
  content: "\e238";
}
.glyphicon-grain:before {
  content: "\e239";
}
.glyphicon-sunglasses:before {
  content: "\e240";
}
.glyphicon-text-size:before {
  content: "\e241";
}
.glyphicon-text-color:before {
  content: "\e242";
}
.glyphicon-text-background:before {
  content: "\e243";
}
.glyphicon-object-align-top:before {
  content: "\e244";
}
.glyphicon-object-align-bottom:before {
  content: "\e245";
}
.glyphicon-object-align-horizontal:before {
  content: "\e246";
}
.glyphicon-object-align-left:before {
  content: "\e247";
}
.glyphicon-object-align-vertical:before {
  content: "\e248";
}
.glyphicon-object-align-right:before {
  content: "\e249";
}
.glyphicon-triangle-right:before {
  content: "\e250";
}
.glyphicon-triangle-left:before {
  content: "\e251";
}
.glyphicon-triangle-bottom:before {
  content: "\e252";
}
.glyphicon-triangle-top:before {
  content: "\e253";
}
.glyphicon-console:before {
  content: "\e254";
}
.glyphicon-superscript:before {
  content: "\e255";
}
.glyphicon-subscript:before {
  content: "\e256";
}
.glyphicon-menu-left:before {
  content: "\e257";
}
.glyphicon-menu-right:before {
  content: "\e258";
}
.glyphicon-menu-down:before {
  content: "\e259";
}
.glyphicon-menu-up:before {
  content: "\e260";
}
* {
  -webkit-box-sizing: border-box;
  -moz-box-sizing: border-box;
  box-sizing: border-box;
}
*:before,
*:after {
  -webkit-box-sizing: border-box;
  -moz-box-sizing: border-box;
  box-sizing: border-box;
}
html {
  font-size: 10px;
  -webkit-tap-highlight-color: rgba(0, 0, 0, 0);
}
body {
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  font-size: 13px;
  line-height: 1.42857143;
  color: #000;
  background-color: #fff;
}
input,
button,
select,
textarea {
  font-family: inherit;
  font-size: inherit;
  line-height: inherit;
}
a {
  color: #337ab7;
  text-decoration: none;
}
a:hover,
a:focus {
  color: #23527c;
  text-decoration: underline;
}
a:focus {
  outline: 5px auto -webkit-focus-ring-color;
  outline-offset: -2px;
}
figure {
  margin: 0;
}
img {
  vertical-align: middle;
}
.img-responsive,
.thumbnail > img,
.thumbnail a > img,
.carousel-inner > .item > img,
.carousel-inner > .item > a > img {
  display: block;
  max-width: 100%;
  height: auto;
}
.img-rounded {
  border-radius: 3px;
}
.img-thumbnail {
  padding: 4px;
  line-height: 1.42857143;
  background-color: #fff;
  border: 1px solid #ddd;
  border-radius: 2px;
  -webkit-transition: all 0.2s ease-in-out;
  -o-transition: all 0.2s ease-in-out;
  transition: all 0.2s ease-in-out;
  display: inline-block;
  max-width: 100%;
  height: auto;
}
.img-circle {
  border-radius: 50%;
}
hr {
  margin-top: 18px;
  margin-bottom: 18px;
  border: 0;
  border-top: 1px solid #eeeeee;
}
.sr-only {
  position: absolute;
  width: 1px;
  height: 1px;
  margin: -1px;
  padding: 0;
  overflow: hidden;
  clip: rect(0, 0, 0, 0);
  border: 0;
}
.sr-only-focusable:active,
.sr-only-focusable:focus {
  position: static;
  width: auto;
  height: auto;
  margin: 0;
  overflow: visible;
  clip: auto;
}
[role="button"] {
  cursor: pointer;
}
h1,
h2,
h3,
h4,
h5,
h6,
.h1,
.h2,
.h3,
.h4,
.h5,
.h6 {
  font-family: inherit;
  font-weight: 500;
  line-height: 1.1;
  color: inherit;
}
h1 small,
h2 small,
h3 small,
h4 small,
h5 small,
h6 small,
.h1 small,
.h2 small,
.h3 small,
.h4 small,
.h5 small,
.h6 small,
h1 .small,
h2 .small,
h3 .small,
h4 .small,
h5 .small,
h6 .small,
.h1 .small,
.h2 .small,
.h3 .small,
.h4 .small,
.h5 .small,
.h6 .small {
  font-weight: normal;
  line-height: 1;
  color: #777777;
}
h1,
.h1,
h2,
.h2,
h3,
.h3 {
  margin-top: 18px;
  margin-bottom: 9px;
}
h1 small,
.h1 small,
h2 small,
.h2 small,
h3 small,
.h3 small,
h1 .small,
.h1 .small,
h2 .small,
.h2 .small,
h3 .small,
.h3 .small {
  font-size: 65%;
}
h4,
.h4,
h5,
.h5,
h6,
.h6 {
  margin-top: 9px;
  margin-bottom: 9px;
}
h4 small,
.h4 small,
h5 small,
.h5 small,
h6 small,
.h6 small,
h4 .small,
.h4 .small,
h5 .small,
.h5 .small,
h6 .small,
.h6 .small {
  font-size: 75%;
}
h1,
.h1 {
  font-size: 33px;
}
h2,
.h2 {
  font-size: 27px;
}
h3,
.h3 {
  font-size: 23px;
}
h4,
.h4 {
  font-size: 17px;
}
h5,
.h5 {
  font-size: 13px;
}
h6,
.h6 {
  font-size: 12px;
}
p {
  margin: 0 0 9px;
}
.lead {
  margin-bottom: 18px;
  font-size: 14px;
  font-weight: 300;
  line-height: 1.4;
}
@media (min-width: 768px) {
  .lead {
    font-size: 19.5px;
  }
}
small,
.small {
  font-size: 92%;
}
mark,
.mark {
  background-color: #fcf8e3;
  padding: .2em;
}
.text-left {
  text-align: left;
}
.text-right {
  text-align: right;
}
.text-center {
  text-align: center;
}
.text-justify {
  text-align: justify;
}
.text-nowrap {
  white-space: nowrap;
}
.text-lowercase {
  text-transform: lowercase;
}
.text-uppercase {
  text-transform: uppercase;
}
.text-capitalize {
  text-transform: capitalize;
}
.text-muted {
  color: #777777;
}
.text-primary {
  color: #337ab7;
}
a.text-primary:hover,
a.text-primary:focus {
  color: #286090;
}
.text-success {
  color: #3c763d;
}
a.text-success:hover,
a.text-success:focus {
  color: #2b542c;
}
.text-info {
  color: #31708f;
}
a.text-info:hover,
a.text-info:focus {
  color: #245269;
}
.text-warning {
  color: #8a6d3b;
}
a.text-warning:hover,
a.text-warning:focus {
  color: #66512c;
}
.text-danger {
  color: #a94442;
}
a.text-danger:hover,
a.text-danger:focus {
  color: #843534;
}
.bg-primary {
  color: #fff;
  background-color: #337ab7;
}
a.bg-primary:hover,
a.bg-primary:focus {
  background-color: #286090;
}
.bg-success {
  background-color: #dff0d8;
}
a.bg-success:hover,
a.bg-success:focus {
  background-color: #c1e2b3;
}
.bg-info {
  background-color: #d9edf7;
}
a.bg-info:hover,
a.bg-info:focus {
  background-color: #afd9ee;
}
.bg-warning {
  background-color: #fcf8e3;
}
a.bg-warning:hover,
a.bg-warning:focus {
  background-color: #f7ecb5;
}
.bg-danger {
  background-color: #f2dede;
}
a.bg-danger:hover,
a.bg-danger:focus {
  background-color: #e4b9b9;
}
.page-header {
  padding-bottom: 8px;
  margin: 36px 0 18px;
  border-bottom: 1px solid #eeeeee;
}
ul,
ol {
  margin-top: 0;
  margin-bottom: 9px;
}
ul ul,
ol ul,
ul ol,
ol ol {
  margin-bottom: 0;
}
.list-unstyled {
  padding-left: 0;
  list-style: none;
}
.list-inline {
  padding-left: 0;
  list-style: none;
  margin-left: -5px;
}
.list-inline > li {
  display: inline-block;
  padding-left: 5px;
  padding-right: 5px;
}
dl {
  margin-top: 0;
  margin-bottom: 18px;
}
dt,
dd {
  line-height: 1.42857143;
}
dt {
  font-weight: bold;
}
dd {
  margin-left: 0;
}
@media (min-width: 541px) {
  .dl-horizontal dt {
    float: left;
    width: 160px;
    clear: left;
    text-align: right;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
  .dl-horizontal dd {
    margin-left: 180px;
  }
}
abbr[title],
abbr[data-original-title] {
  cursor: help;
  border-bottom: 1px dotted #777777;
}
.initialism {
  font-size: 90%;
  text-transform: uppercase;
}
blockquote {
  padding: 9px 18px;
  margin: 0 0 18px;
  font-size: inherit;
  border-left: 5px solid #eeeeee;
}
blockquote p:last-child,
blockquote ul:last-child,
blockquote ol:last-child {
  margin-bottom: 0;
}
blockquote footer,
blockquote small,
blockquote .small {
  display: block;
  font-size: 80%;
  line-height: 1.42857143;
  color: #777777;
}
blockquote footer:before,
blockquote small:before,
blockquote .small:before {
  content: '\2014 \00A0';
}
.blockquote-reverse,
blockquote.pull-right {
  padding-right: 15px;
  padding-left: 0;
  border-right: 5px solid #eeeeee;
  border-left: 0;
  text-align: right;
}
.blockquote-reverse footer:before,
blockquote.pull-right footer:before,
.blockquote-reverse small:before,
blockquote.pull-right small:before,
.blockquote-reverse .small:before,
blockquote.pull-right .small:before {
  content: '';
}
.blockquote-reverse footer:after,
blockquote.pull-right footer:after,
.blockquote-reverse small:after,
blockquote.pull-right small:after,
.blockquote-reverse .small:after,
blockquote.pull-right .small:after {
  content: '\00A0 \2014';
}
address {
  margin-bottom: 18px;
  font-style: normal;
  line-height: 1.42857143;
}
code,
kbd,
pre,
samp {
  font-family: monospace;
}
code {
  padding: 2px 4px;
  font-size: 90%;
  color: #c7254e;
  background-color: #f9f2f4;
  border-radius: 2px;
}
kbd {
  padding: 2px 4px;
  font-size: 90%;
  color: #888;
  background-color: transparent;
  border-radius: 1px;
  box-shadow: inset 0 -1px 0 rgba(0, 0, 0, 0.25);
}
kbd kbd {
  padding: 0;
  font-size: 100%;
  font-weight: bold;
  box-shadow: none;
}
pre {
  display: block;
  padding: 8.5px;
  margin: 0 0 9px;
  font-size: 12px;
  line-height: 1.42857143;
  word-break: break-all;
  word-wrap: break-word;
  color: #333333;
  background-color: #f5f5f5;
  border: 1px solid #ccc;
  border-radius: 2px;
}
pre code {
  padding: 0;
  font-size: inherit;
  color: inherit;
  white-space: pre-wrap;
  background-color: transparent;
  border-radius: 0;
}
.pre-scrollable {
  max-height: 340px;
  overflow-y: scroll;
}
.container {
  margin-right: auto;
  margin-left: auto;
  padding-left: 0px;
  padding-right: 0px;
}
@media (min-width: 768px) {
  .container {
    width: 768px;
  }
}
@media (min-width: 992px) {
  .container {
    width: 940px;
  }
}
@media (min-width: 1200px) {
  .container {
    width: 1140px;
  }
}
.container-fluid {
  margin-right: auto;
  margin-left: auto;
  padding-left: 0px;
  padding-right: 0px;
}
.row {
  margin-left: 0px;
  margin-right: 0px;
}
.col-xs-1, .col-sm-1, .col-md-1, .col-lg-1, .col-xs-2, .col-sm-2, .col-md-2, .col-lg-2, .col-xs-3, .col-sm-3, .col-md-3, .col-lg-3, .col-xs-4, .col-sm-4, .col-md-4, .col-lg-4, .col-xs-5, .col-sm-5, .col-md-5, .col-lg-5, .col-xs-6, .col-sm-6, .col-md-6, .col-lg-6, .col-xs-7, .col-sm-7, .col-md-7, .col-lg-7, .col-xs-8, .col-sm-8, .col-md-8, .col-lg-8, .col-xs-9, .col-sm-9, .col-md-9, .col-lg-9, .col-xs-10, .col-sm-10, .col-md-10, .col-lg-10, .col-xs-11, .col-sm-11, .col-md-11, .col-lg-11, .col-xs-12, .col-sm-12, .col-md-12, .col-lg-12 {
  position: relative;
  min-height: 1px;
  padding-left: 0px;
  padding-right: 0px;
}
.col-xs-1, .col-xs-2, .col-xs-3, .col-xs-4, .col-xs-5, .col-xs-6, .col-xs-7, .col-xs-8, .col-xs-9, .col-xs-10, .col-xs-11, .col-xs-12 {
  float: left;
}
.col-xs-12 {
  width: 100%;
}
.col-xs-11 {
  width: 91.66666667%;
}
.col-xs-10 {
  width: 83.33333333%;
}
.col-xs-9 {
  width: 75%;
}
.col-xs-8 {
  width: 66.66666667%;
}
.col-xs-7 {
  width: 58.33333333%;
}
.col-xs-6 {
  width: 50%;
}
.col-xs-5 {
  width: 41.66666667%;
}
.col-xs-4 {
  width: 33.33333333%;
}
.col-xs-3 {
  width: 25%;
}
.col-xs-2 {
  width: 16.66666667%;
}
.col-xs-1 {
  width: 8.33333333%;
}
.col-xs-pull-12 {
  right: 100%;
}
.col-xs-pull-11 {
  right: 91.66666667%;
}
.col-xs-pull-10 {
  right: 83.33333333%;
}
.col-xs-pull-9 {
  right: 75%;
}
.col-xs-pull-8 {
  right: 66.66666667%;
}
.col-xs-pull-7 {
  right: 58.33333333%;
}
.col-xs-pull-6 {
  right: 50%;
}
.col-xs-pull-5 {
  right: 41.66666667%;
}
.col-xs-pull-4 {
  right: 33.33333333%;
}
.col-xs-pull-3 {
  right: 25%;
}
.col-xs-pull-2 {
  right: 16.66666667%;
}
.col-xs-pull-1 {
  right: 8.33333333%;
}
.col-xs-pull-0 {
  right: auto;
}
.col-xs-push-12 {
  left: 100%;
}
.col-xs-push-11 {
  left: 91.66666667%;
}
.col-xs-push-10 {
  left: 83.33333333%;
}
.col-xs-push-9 {
  left: 75%;
}
.col-xs-push-8 {
  left: 66.66666667%;
}
.col-xs-push-7 {
  left: 58.33333333%;
}
.col-xs-push-6 {
  left: 50%;
}
.col-xs-push-5 {
  left: 41.66666667%;
}
.col-xs-push-4 {
  left: 33.33333333%;
}
.col-xs-push-3 {
  left: 25%;
}
.col-xs-push-2 {
  left: 16.66666667%;
}
.col-xs-push-1 {
  left: 8.33333333%;
}
.col-xs-push-0 {
  left: auto;
}
.col-xs-offset-12 {
  margin-left: 100%;
}
.col-xs-offset-11 {
  margin-left: 91.66666667%;
}
.col-xs-offset-10 {
  margin-left: 83.33333333%;
}
.col-xs-offset-9 {
  margin-left: 75%;
}
.col-xs-offset-8 {
  margin-left: 66.66666667%;
}
.col-xs-offset-7 {
  margin-left: 58.33333333%;
}
.col-xs-offset-6 {
  margin-left: 50%;
}
.col-xs-offset-5 {
  margin-left: 41.66666667%;
}
.col-xs-offset-4 {
  margin-left: 33.33333333%;
}
.col-xs-offset-3 {
  margin-left: 25%;
}
.col-xs-offset-2 {
  margin-left: 16.66666667%;
}
.col-xs-offset-1 {
  margin-left: 8.33333333%;
}
.col-xs-offset-0 {
  margin-left: 0%;
}
@media (min-width: 768px) {
  .col-sm-1, .col-sm-2, .col-sm-3, .col-sm-4, .col-sm-5, .col-sm-6, .col-sm-7, .col-sm-8, .col-sm-9, .col-sm-10, .col-sm-11, .col-sm-12 {
    float: left;
  }
  .col-sm-12 {
    width: 100%;
  }
  .col-sm-11 {
    width: 91.66666667%;
  }
  .col-sm-10 {
    width: 83.33333333%;
  }
  .col-sm-9 {
    width: 75%;
  }
  .col-sm-8 {
    width: 66.66666667%;
  }
  .col-sm-7 {
    width: 58.33333333%;
  }
  .col-sm-6 {
    width: 50%;
  }
  .col-sm-5 {
    width: 41.66666667%;
  }
  .col-sm-4 {
    width: 33.33333333%;
  }
  .col-sm-3 {
    width: 25%;
  }
  .col-sm-2 {
    width: 16.66666667%;
  }
  .col-sm-1 {
    width: 8.33333333%;
  }
  .col-sm-pull-12 {
    right: 100%;
  }
  .col-sm-pull-11 {
    right: 91.66666667%;
  }
  .col-sm-pull-10 {
    right: 83.33333333%;
  }
  .col-sm-pull-9 {
    right: 75%;
  }
  .col-sm-pull-8 {
    right: 66.66666667%;
  }
  .col-sm-pull-7 {
    right: 58.33333333%;
  }
  .col-sm-pull-6 {
    right: 50%;
  }
  .col-sm-pull-5 {
    right: 41.66666667%;
  }
  .col-sm-pull-4 {
    right: 33.33333333%;
  }
  .col-sm-pull-3 {
    right: 25%;
  }
  .col-sm-pull-2 {
    right: 16.66666667%;
  }
  .col-sm-pull-1 {
    right: 8.33333333%;
  }
  .col-sm-pull-0 {
    right: auto;
  }
  .col-sm-push-12 {
    left: 100%;
  }
  .col-sm-push-11 {
    left: 91.66666667%;
  }
  .col-sm-push-10 {
    left: 83.33333333%;
  }
  .col-sm-push-9 {
    left: 75%;
  }
  .col-sm-push-8 {
    left: 66.66666667%;
  }
  .col-sm-push-7 {
    left: 58.33333333%;
  }
  .col-sm-push-6 {
    left: 50%;
  }
  .col-sm-push-5 {
    left: 41.66666667%;
  }
  .col-sm-push-4 {
    left: 33.33333333%;
  }
  .col-sm-push-3 {
    left: 25%;
  }
  .col-sm-push-2 {
    left: 16.66666667%;
  }
  .col-sm-push-1 {
    left: 8.33333333%;
  }
  .col-sm-push-0 {
    left: auto;
  }
  .col-sm-offset-12 {
    margin-left: 100%;
  }
  .col-sm-offset-11 {
    margin-left: 91.66666667%;
  }
  .col-sm-offset-10 {
    margin-left: 83.33333333%;
  }
  .col-sm-offset-9 {
    margin-left: 75%;
  }
  .col-sm-offset-8 {
    margin-left: 66.66666667%;
  }
  .col-sm-offset-7 {
    margin-left: 58.33333333%;
  }
  .col-sm-offset-6 {
    margin-left: 50%;
  }
  .col-sm-offset-5 {
    margin-left: 41.66666667%;
  }
  .col-sm-offset-4 {
    margin-left: 33.33333333%;
  }
  .col-sm-offset-3 {
    margin-left: 25%;
  }
  .col-sm-offset-2 {
    margin-left: 16.66666667%;
  }
  .col-sm-offset-1 {
    margin-left: 8.33333333%;
  }
  .col-sm-offset-0 {
    margin-left: 0%;
  }
}
@media (min-width: 992px) {
  .col-md-1, .col-md-2, .col-md-3, .col-md-4, .col-md-5, .col-md-6, .col-md-7, .col-md-8, .col-md-9, .col-md-10, .col-md-11, .col-md-12 {
    float: left;
  }
  .col-md-12 {
    width: 100%;
  }
  .col-md-11 {
    width: 91.66666667%;
  }
  .col-md-10 {
    width: 83.33333333%;
  }
  .col-md-9 {
    width: 75%;
  }
  .col-md-8 {
    width: 66.66666667%;
  }
  .col-md-7 {
    width: 58.33333333%;
  }
  .col-md-6 {
    width: 50%;
  }
  .col-md-5 {
    width: 41.66666667%;
  }
  .col-md-4 {
    width: 33.33333333%;
  }
  .col-md-3 {
    width: 25%;
  }
  .col-md-2 {
    width: 16.66666667%;
  }
  .col-md-1 {
    width: 8.33333333%;
  }
  .col-md-pull-12 {
    right: 100%;
  }
  .col-md-pull-11 {
    right: 91.66666667%;
  }
  .col-md-pull-10 {
    right: 83.33333333%;
  }
  .col-md-pull-9 {
    right: 75%;
  }
  .col-md-pull-8 {
    right: 66.66666667%;
  }
  .col-md-pull-7 {
    right: 58.33333333%;
  }
  .col-md-pull-6 {
    right: 50%;
  }
  .col-md-pull-5 {
    right: 41.66666667%;
  }
  .col-md-pull-4 {
    right: 33.33333333%;
  }
  .col-md-pull-3 {
    right: 25%;
  }
  .col-md-pull-2 {
    right: 16.66666667%;
  }
  .col-md-pull-1 {
    right: 8.33333333%;
  }
  .col-md-pull-0 {
    right: auto;
  }
  .col-md-push-12 {
    left: 100%;
  }
  .col-md-push-11 {
    left: 91.66666667%;
  }
  .col-md-push-10 {
    left: 83.33333333%;
  }
  .col-md-push-9 {
    left: 75%;
  }
  .col-md-push-8 {
    left: 66.66666667%;
  }
  .col-md-push-7 {
    left: 58.33333333%;
  }
  .col-md-push-6 {
    left: 50%;
  }
  .col-md-push-5 {
    left: 41.66666667%;
  }
  .col-md-push-4 {
    left: 33.33333333%;
  }
  .col-md-push-3 {
    left: 25%;
  }
  .col-md-push-2 {
    left: 16.66666667%;
  }
  .col-md-push-1 {
    left: 8.33333333%;
  }
  .col-md-push-0 {
    left: auto;
  }
  .col-md-offset-12 {
    margin-left: 100%;
  }
  .col-md-offset-11 {
    margin-left: 91.66666667%;
  }
  .col-md-offset-10 {
    margin-left: 83.33333333%;
  }
  .col-md-offset-9 {
    margin-left: 75%;
  }
  .col-md-offset-8 {
    margin-left: 66.66666667%;
  }
  .col-md-offset-7 {
    margin-left: 58.33333333%;
  }
  .col-md-offset-6 {
    margin-left: 50%;
  }
  .col-md-offset-5 {
    margin-left: 41.66666667%;
  }
  .col-md-offset-4 {
    margin-left: 33.33333333%;
  }
  .col-md-offset-3 {
    margin-left: 25%;
  }
  .col-md-offset-2 {
    margin-left: 16.66666667%;
  }
  .col-md-offset-1 {
    margin-left: 8.33333333%;
  }
  .col-md-offset-0 {
    margin-left: 0%;
  }
}
@media (min-width: 1200px) {
  .col-lg-1, .col-lg-2, .col-lg-3, .col-lg-4, .col-lg-5, .col-lg-6, .col-lg-7, .col-lg-8, .col-lg-9, .col-lg-10, .col-lg-11, .col-lg-12 {
    float: left;
  }
  .col-lg-12 {
    width: 100%;
  }
  .col-lg-11 {
    width: 91.66666667%;
  }
  .col-lg-10 {
    width: 83.33333333%;
  }
  .col-lg-9 {
    width: 75%;
  }
  .col-lg-8 {
    width: 66.66666667%;
  }
  .col-lg-7 {
    width: 58.33333333%;
  }
  .col-lg-6 {
    width: 50%;
  }
  .col-lg-5 {
    width: 41.66666667%;
  }
  .col-lg-4 {
    width: 33.33333333%;
  }
  .col-lg-3 {
    width: 25%;
  }
  .col-lg-2 {
    width: 16.66666667%;
  }
  .col-lg-1 {
    width: 8.33333333%;
  }
  .col-lg-pull-12 {
    right: 100%;
  }
  .col-lg-pull-11 {
    right: 91.66666667%;
  }
  .col-lg-pull-10 {
    right: 83.33333333%;
  }
  .col-lg-pull-9 {
    right: 75%;
  }
  .col-lg-pull-8 {
    right: 66.66666667%;
  }
  .col-lg-pull-7 {
    right: 58.33333333%;
  }
  .col-lg-pull-6 {
    right: 50%;
  }
  .col-lg-pull-5 {
    right: 41.66666667%;
  }
  .col-lg-pull-4 {
    right: 33.33333333%;
  }
  .col-lg-pull-3 {
    right: 25%;
  }
  .col-lg-pull-2 {
    right: 16.66666667%;
  }
  .col-lg-pull-1 {
    right: 8.33333333%;
  }
  .col-lg-pull-0 {
    right: auto;
  }
  .col-lg-push-12 {
    left: 100%;
  }
  .col-lg-push-11 {
    left: 91.66666667%;
  }
  .col-lg-push-10 {
    left: 83.33333333%;
  }
  .col-lg-push-9 {
    left: 75%;
  }
  .col-lg-push-8 {
    left: 66.66666667%;
  }
  .col-lg-push-7 {
    left: 58.33333333%;
  }
  .col-lg-push-6 {
    left: 50%;
  }
  .col-lg-push-5 {
    left: 41.66666667%;
  }
  .col-lg-push-4 {
    left: 33.33333333%;
  }
  .col-lg-push-3 {
    left: 25%;
  }
  .col-lg-push-2 {
    left: 16.66666667%;
  }
  .col-lg-push-1 {
    left: 8.33333333%;
  }
  .col-lg-push-0 {
    left: auto;
  }
  .col-lg-offset-12 {
    margin-left: 100%;
  }
  .col-lg-offset-11 {
    margin-left: 91.66666667%;
  }
  .col-lg-offset-10 {
    margin-left: 83.33333333%;
  }
  .col-lg-offset-9 {
    margin-left: 75%;
  }
  .col-lg-offset-8 {
    margin-left: 66.66666667%;
  }
  .col-lg-offset-7 {
    margin-left: 58.33333333%;
  }
  .col-lg-offset-6 {
    margin-left: 50%;
  }
  .col-lg-offset-5 {
    margin-left: 41.66666667%;
  }
  .col-lg-offset-4 {
    margin-left: 33.33333333%;
  }
  .col-lg-offset-3 {
    margin-left: 25%;
  }
  .col-lg-offset-2 {
    margin-left: 16.66666667%;
  }
  .col-lg-offset-1 {
    margin-left: 8.33333333%;
  }
  .col-lg-offset-0 {
    margin-left: 0%;
  }
}
table {
  background-color: transparent;
}
caption {
  padding-top: 8px;
  padding-bottom: 8px;
  color: #777777;
  text-align: left;
}
th {
  text-align: left;
}
.table {
  width: 100%;
  max-width: 100%;
  margin-bottom: 18px;
}
.table > thead > tr > th,
.table > tbody > tr > th,
.table > tfoot > tr > th,
.table > thead > tr > td,
.table > tbody > tr > td,
.table > tfoot > tr > td {
  padding: 8px;
  line-height: 1.42857143;
  vertical-align: top;
  border-top: 1px solid #ddd;
}
.table > thead > tr > th {
  vertical-align: bottom;
  border-bottom: 2px solid #ddd;
}
.table > caption + thead > tr:first-child > th,
.table > colgroup + thead > tr:first-child > th,
.table > thead:first-child > tr:first-child > th,
.table > caption + thead > tr:first-child > td,
.table > colgroup + thead > tr:first-child > td,
.table > thead:first-child > tr:first-child > td {
  border-top: 0;
}
.table > tbody + tbody {
  border-top: 2px solid #ddd;
}
.table .table {
  background-color: #fff;
}
.table-condensed > thead > tr > th,
.table-condensed > tbody > tr > th,
.table-condensed > tfoot > tr > th,
.table-condensed > thead > tr > td,
.table-condensed > tbody > tr > td,
.table-condensed > tfoot > tr > td {
  padding: 5px;
}
.table-bordered {
  border: 1px solid #ddd;
}
.table-bordered > thead > tr > th,
.table-bordered > tbody > tr > th,
.table-bordered > tfoot > tr > th,
.table-bordered > thead > tr > td,
.table-bordered > tbody > tr > td,
.table-bordered > tfoot > tr > td {
  border: 1px solid #ddd;
}
.table-bordered > thead > tr > th,
.table-bordered > thead > tr > td {
  border-bottom-width: 2px;
}
.table-striped > tbody > tr:nth-of-type(odd) {
  background-color: #f9f9f9;
}
.table-hover > tbody > tr:hover {
  background-color: #f5f5f5;
}
table col[class*="col-"] {
  position: static;
  float: none;
  display: table-column;
}
table td[class*="col-"],
table th[class*="col-"] {
  position: static;
  float: none;
  display: table-cell;
}
.table > thead > tr > td.active,
.table > tbody > tr > td.active,
.table > tfoot > tr > td.active,
.table > thead > tr > th.active,
.table > tbody > tr > th.active,
.table > tfoot > tr > th.active,
.table > thead > tr.active > td,
.table > tbody > tr.active > td,
.table > tfoot > tr.active > td,
.table > thead > tr.active > th,
.table > tbody > tr.active > th,
.table > tfoot > tr.active > th {
  background-color: #f5f5f5;
}
.table-hover > tbody > tr > td.active:hover,
.table-hover > tbody > tr > th.active:hover,
.table-hover > tbody > tr.active:hover > td,
.table-hover > tbody > tr:hover > .active,
.table-hover > tbody > tr.active:hover > th {
  background-color: #e8e8e8;
}
.table > thead > tr > td.success,
.table > tbody > tr > td.success,
.table > tfoot > tr > td.success,
.table > thead > tr > th.success,
.table > tbody > tr > th.success,
.table > tfoot > tr > th.success,
.table > thead > tr.success > td,
.table > tbody > tr.success > td,
.table > tfoot > tr.success > td,
.table > thead > tr.success > th,
.table > tbody > tr.success > th,
.table > tfoot > tr.success > th {
  background-color: #dff0d8;
}
.table-hover > tbody > tr > td.success:hover,
.table-hover > tbody > tr > th.success:hover,
.table-hover > tbody > tr.success:hover > td,
.table-hover > tbody > tr:hover > .success,
.table-hover > tbody > tr.success:hover > th {
  background-color: #d0e9c6;
}
.table > thead > tr > td.info,
.table > tbody > tr > td.info,
.table > tfoot > tr > td.info,
.table > thead > tr > th.info,
.table > tbody > tr > th.info,
.table > tfoot > tr > th.info,
.table > thead > tr.info > td,
.table > tbody > tr.info > td,
.table > tfoot > tr.info > td,
.table > thead > tr.info > th,
.table > tbody > tr.info > th,
.table > tfoot > tr.info > th {
  background-color: #d9edf7;
}
.table-hover > tbody > tr > td.info:hover,
.table-hover > tbody > tr > th.info:hover,
.table-hover > tbody > tr.info:hover > td,
.table-hover > tbody > tr:hover > .info,
.table-hover > tbody > tr.info:hover > th {
  background-color: #c4e3f3;
}
.table > thead > tr > td.warning,
.table > tbody > tr > td.warning,
.table > tfoot > tr > td.warning,
.table > thead > tr > th.warning,
.table > tbody > tr > th.warning,
.table > tfoot > tr > th.warning,
.table > thead > tr.warning > td,
.table > tbody > tr.warning > td,
.table > tfoot > tr.warning > td,
.table > thead > tr.warning > th,
.table > tbody > tr.warning > th,
.table > tfoot > tr.warning > th {
  background-color: #fcf8e3;
}
.table-hover > tbody > tr > td.warning:hover,
.table-hover > tbody > tr > th.warning:hover,
.table-hover > tbody > tr.warning:hover > td,
.table-hover > tbody > tr:hover > .warning,
.table-hover > tbody > tr.warning:hover > th {
  background-color: #faf2cc;
}
.table > thead > tr > td.danger,
.table > tbody > tr > td.danger,
.table > tfoot > tr > td.danger,
.table > thead > tr > th.danger,
.table > tbody > tr > th.danger,
.table > tfoot > tr > th.danger,
.table > thead > tr.danger > td,
.table > tbody > tr.danger > td,
.table > tfoot > tr.danger > td,
.table > thead > tr.danger > th,
.table > tbody > tr.danger > th,
.table > tfoot > tr.danger > th {
  background-color: #f2dede;
}
.table-hover > tbody > tr > td.danger:hover,
.table-hover > tbody > tr > th.danger:hover,
.table-hover > tbody > tr.danger:hover > td,
.table-hover > tbody > tr:hover > .danger,
.table-hover > tbody > tr.danger:hover > th {
  background-color: #ebcccc;
}
.table-responsive {
  overflow-x: auto;
  min-height: 0.01%;
}
@media screen and (max-width: 767px) {
  .table-responsive {
    width: 100%;
    margin-bottom: 13.5px;
    overflow-y: hidden;
    -ms-overflow-style: -ms-autohiding-scrollbar;
    border: 1px solid #ddd;
  }
  .table-responsive > .table {
    margin-bottom: 0;
  }
  .table-responsive > .table > thead > tr > th,
  .table-responsive > .table > tbody > tr > th,
  .table-responsive > .table > tfoot > tr > th,
  .table-responsive > .table > thead > tr > td,
  .table-responsive > .table > tbody > tr > td,
  .table-responsive > .table > tfoot > tr > td {
    white-space: nowrap;
  }
  .table-responsive > .table-bordered {
    border: 0;
  }
  .table-responsive > .table-bordered > thead > tr > th:first-child,
  .table-responsive > .table-bordered > tbody > tr > th:first-child,
  .table-responsive > .table-bordered > tfoot > tr > th:first-child,
  .table-responsive > .table-bordered > thead > tr > td:first-child,
  .table-responsive > .table-bordered > tbody > tr > td:first-child,
  .table-responsive > .table-bordered > tfoot > tr > td:first-child {
    border-left: 0;
  }
  .table-responsive > .table-bordered > thead > tr > th:last-child,
  .table-responsive > .table-bordered > tbody > tr > th:last-child,
  .table-responsive > .table-bordered > tfoot > tr > th:last-child,
  .table-responsive > .table-bordered > thead > tr > td:last-child,
  .table-responsive > .table-bordered > tbody > tr > td:last-child,
  .table-responsive > .table-bordered > tfoot > tr > td:last-child {
    border-right: 0;
  }
  .table-responsive > .table-bordered > tbody > tr:last-child > th,
  .table-responsive > .table-bordered > tfoot > tr:last-child > th,
  .table-responsive > .table-bordered > tbody > tr:last-child > td,
  .table-responsive > .table-bordered > tfoot > tr:last-child > td {
    border-bottom: 0;
  }
}
fieldset {
  padding: 0;
  margin: 0;
  border: 0;
  min-width: 0;
}
legend {
  display: block;
  width: 100%;
  padding: 0;
  margin-bottom: 18px;
  font-size: 19.5px;
  line-height: inherit;
  color: #333333;
  border: 0;
  border-bottom: 1px solid #e5e5e5;
}
label {
  display: inline-block;
  max-width: 100%;
  margin-bottom: 5px;
  font-weight: bold;
}
input[type="search"] {
  -webkit-box-sizing: border-box;
  -moz-box-sizing: border-box;
  box-sizing: border-box;
}
input[type="radio"],
input[type="checkbox"] {
  margin: 4px 0 0;
  margin-top: 1px \9;
  line-height: normal;
}
input[type="file"] {
  display: block;
}
input[type="range"] {
  display: block;
  width: 100%;
}
select[multiple],
select[size] {
  height: auto;
}
input[type="file"]:focus,
input[type="radio"]:focus,
input[type="checkbox"]:focus {
  outline: 5px auto -webkit-focus-ring-color;
  outline-offset: -2px;
}
output {
  display: block;
  padding-top: 7px;
  font-size: 13px;
  line-height: 1.42857143;
  color: #555555;
}
.form-control {
  display: block;
  width: 100%;
  height: 32px;
  padding: 6px 12px;
  font-size: 13px;
  line-height: 1.42857143;
  color: #555555;
  background-color: #fff;
  background-image: none;
  border: 1px solid #ccc;
  border-radius: 2px;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  -webkit-transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
  -o-transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
  transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
}
.form-control:focus {
  border-color: #66afe9;
  outline: 0;
  -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.075), 0 0 8px rgba(102, 175, 233, 0.6);
  box-shadow: inset 0 1px 1px rgba(0,0,0,.075), 0 0 8px rgba(102, 175, 233, 0.6);
}
.form-control::-moz-placeholder {
  color: #999;
  opacity: 1;
}
.form-control:-ms-input-placeholder {
  color: #999;
}
.form-control::-webkit-input-placeholder {
  color: #999;
}
.form-control::-ms-expand {
  border: 0;
  background-color: transparent;
}
.form-control[disabled],
.form-control[readonly],
fieldset[disabled] .form-control {
  background-color: #eeeeee;
  opacity: 1;
}
.form-control[disabled],
fieldset[disabled] .form-control {
  cursor: not-allowed;
}
textarea.form-control {
  height: auto;
}
input[type="search"] {
  -webkit-appearance: none;
}
@media screen and (-webkit-min-device-pixel-ratio: 0) {
  input[type="date"].form-control,
  input[type="time"].form-control,
  input[type="datetime-local"].form-control,
  input[type="month"].form-control {
    line-height: 32px;
  }
  input[type="date"].input-sm,
  input[type="time"].input-sm,
  input[type="datetime-local"].input-sm,
  input[type="month"].input-sm,
  .input-group-sm input[type="date"],
  .input-group-sm input[type="time"],
  .input-group-sm input[type="datetime-local"],
  .input-group-sm input[type="month"] {
    line-height: 30px;
  }
  input[type="date"].input-lg,
  input[type="time"].input-lg,
  input[type="datetime-local"].input-lg,
  input[type="month"].input-lg,
  .input-group-lg input[type="date"],
  .input-group-lg input[type="time"],
  .input-group-lg input[type="datetime-local"],
  .input-group-lg input[type="month"] {
    line-height: 45px;
  }
}
.form-group {
  margin-bottom: 15px;
}
.radio,
.checkbox {
  position: relative;
  display: block;
  margin-top: 10px;
  margin-bottom: 10px;
}
.radio label,
.checkbox label {
  min-height: 18px;
  padding-left: 20px;
  margin-bottom: 0;
  font-weight: normal;
  cursor: pointer;
}
.radio input[type="radio"],
.radio-inline input[type="radio"],
.checkbox input[type="checkbox"],
.checkbox-inline input[type="checkbox"] {
  position: absolute;
  margin-left: -20px;
  margin-top: 4px \9;
}
.radio + .radio,
.checkbox + .checkbox {
  margin-top: -5px;
}
.radio-inline,
.checkbox-inline {
  position: relative;
  display: inline-block;
  padding-left: 20px;
  margin-bottom: 0;
  vertical-align: middle;
  font-weight: normal;
  cursor: pointer;
}
.radio-inline + .radio-inline,
.checkbox-inline + .checkbox-inline {
  margin-top: 0;
  margin-left: 10px;
}
input[type="radio"][disabled],
input[type="checkbox"][disabled],
input[type="radio"].disabled,
input[type="checkbox"].disabled,
fieldset[disabled] input[type="radio"],
fieldset[disabled] input[type="checkbox"] {
  cursor: not-allowed;
}
.radio-inline.disabled,
.checkbox-inline.disabled,
fieldset[disabled] .radio-inline,
fieldset[disabled] .checkbox-inline {
  cursor: not-allowed;
}
.radio.disabled label,
.checkbox.disabled label,
fieldset[disabled] .radio label,
fieldset[disabled] .checkbox label {
  cursor: not-allowed;
}
.form-control-static {
  padding-top: 7px;
  padding-bottom: 7px;
  margin-bottom: 0;
  min-height: 31px;
}
.form-control-static.input-lg,
.form-control-static.input-sm {
  padding-left: 0;
  padding-right: 0;
}
.input-sm {
  height: 30px;
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
}
select.input-sm {
  height: 30px;
  line-height: 30px;
}
textarea.input-sm,
select[multiple].input-sm {
  height: auto;
}
.form-group-sm .form-control {
  height: 30px;
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
}
.form-group-sm select.form-control {
  height: 30px;
  line-height: 30px;
}
.form-group-sm textarea.form-control,
.form-group-sm select[multiple].form-control {
  height: auto;
}
.form-group-sm .form-control-static {
  height: 30px;
  min-height: 30px;
  padding: 6px 10px;
  font-size: 12px;
  line-height: 1.5;
}
.input-lg {
  height: 45px;
  padding: 10px 16px;
  font-size: 17px;
  line-height: 1.3333333;
  border-radius: 3px;
}
select.input-lg {
  height: 45px;
  line-height: 45px;
}
textarea.input-lg,
select[multiple].input-lg {
  height: auto;
}
.form-group-lg .form-control {
  height: 45px;
  padding: 10px 16px;
  font-size: 17px;
  line-height: 1.3333333;
  border-radius: 3px;
}
.form-group-lg select.form-control {
  height: 45px;
  line-height: 45px;
}
.form-group-lg textarea.form-control,
.form-group-lg select[multiple].form-control {
  height: auto;
}
.form-group-lg .form-control-static {
  height: 45px;
  min-height: 35px;
  padding: 11px 16px;
  font-size: 17px;
  line-height: 1.3333333;
}
.has-feedback {
  position: relative;
}
.has-feedback .form-control {
  padding-right: 40px;
}
.form-control-feedback {
  position: absolute;
  top: 0;
  right: 0;
  z-index: 2;
  display: block;
  width: 32px;
  height: 32px;
  line-height: 32px;
  text-align: center;
  pointer-events: none;
}
.input-lg + .form-control-feedback,
.input-group-lg + .form-control-feedback,
.form-group-lg .form-control + .form-control-feedback {
  width: 45px;
  height: 45px;
  line-height: 45px;
}
.input-sm + .form-control-feedback,
.input-group-sm + .form-control-feedback,
.form-group-sm .form-control + .form-control-feedback {
  width: 30px;
  height: 30px;
  line-height: 30px;
}
.has-success .help-block,
.has-success .control-label,
.has-success .radio,
.has-success .checkbox,
.has-success .radio-inline,
.has-success .checkbox-inline,
.has-success.radio label,
.has-success.checkbox label,
.has-success.radio-inline label,
.has-success.checkbox-inline label {
  color: #3c763d;
}
.has-success .form-control {
  border-color: #3c763d;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
}
.has-success .form-control:focus {
  border-color: #2b542c;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #67b168;
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #67b168;
}
.has-success .input-group-addon {
  color: #3c763d;
  border-color: #3c763d;
  background-color: #dff0d8;
}
.has-success .form-control-feedback {
  color: #3c763d;
}
.has-warning .help-block,
.has-warning .control-label,
.has-warning .radio,
.has-warning .checkbox,
.has-warning .radio-inline,
.has-warning .checkbox-inline,
.has-warning.radio label,
.has-warning.checkbox label,
.has-warning.radio-inline label,
.has-warning.checkbox-inline label {
  color: #8a6d3b;
}
.has-warning .form-control {
  border-color: #8a6d3b;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
}
.has-warning .form-control:focus {
  border-color: #66512c;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #c0a16b;
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #c0a16b;
}
.has-warning .input-group-addon {
  color: #8a6d3b;
  border-color: #8a6d3b;
  background-color: #fcf8e3;
}
.has-warning .form-control-feedback {
  color: #8a6d3b;
}
.has-error .help-block,
.has-error .control-label,
.has-error .radio,
.has-error .checkbox,
.has-error .radio-inline,
.has-error .checkbox-inline,
.has-error.radio label,
.has-error.checkbox label,
.has-error.radio-inline label,
.has-error.checkbox-inline label {
  color: #a94442;
}
.has-error .form-control {
  border-color: #a94442;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
}
.has-error .form-control:focus {
  border-color: #843534;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #ce8483;
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075), 0 0 6px #ce8483;
}
.has-error .input-group-addon {
  color: #a94442;
  border-color: #a94442;
  background-color: #f2dede;
}
.has-error .form-control-feedback {
  color: #a94442;
}
.has-feedback label ~ .form-control-feedback {
  top: 23px;
}
.has-feedback label.sr-only ~ .form-control-feedback {
  top: 0;
}
.help-block {
  display: block;
  margin-top: 5px;
  margin-bottom: 10px;
  color: #404040;
}
@media (min-width: 768px) {
  .form-inline .form-group {
    display: inline-block;
    margin-bottom: 0;
    vertical-align: middle;
  }
  .form-inline .form-control {
    display: inline-block;
    width: auto;
    vertical-align: middle;
  }
  .form-inline .form-control-static {
    display: inline-block;
  }
  .form-inline .input-group {
    display: inline-table;
    vertical-align: middle;
  }
  .form-inline .input-group .input-group-addon,
  .form-inline .input-group .input-group-btn,
  .form-inline .input-group .form-control {
    width: auto;
  }
  .form-inline .input-group > .form-control {
    width: 100%;
  }
  .form-inline .control-label {
    margin-bottom: 0;
    vertical-align: middle;
  }
  .form-inline .radio,
  .form-inline .checkbox {
    display: inline-block;
    margin-top: 0;
    margin-bottom: 0;
    vertical-align: middle;
  }
  .form-inline .radio label,
  .form-inline .checkbox label {
    padding-left: 0;
  }
  .form-inline .radio input[type="radio"],
  .form-inline .checkbox input[type="checkbox"] {
    position: relative;
    margin-left: 0;
  }
  .form-inline .has-feedback .form-control-feedback {
    top: 0;
  }
}
.form-horizontal .radio,
.form-horizontal .checkbox,
.form-horizontal .radio-inline,
.form-horizontal .checkbox-inline {
  margin-top: 0;
  margin-bottom: 0;
  padding-top: 7px;
}
.form-horizontal .radio,
.form-horizontal .checkbox {
  min-height: 25px;
}
.form-horizontal .form-group {
  margin-left: 0px;
  margin-right: 0px;
}
@media (min-width: 768px) {
  .form-horizontal .control-label {
    text-align: right;
    margin-bottom: 0;
    padding-top: 7px;
  }
}
.form-horizontal .has-feedback .form-control-feedback {
  right: 0px;
}
@media (min-width: 768px) {
  .form-horizontal .form-group-lg .control-label {
    padding-top: 11px;
    font-size: 17px;
  }
}
@media (min-width: 768px) {
  .form-horizontal .form-group-sm .control-label {
    padding-top: 6px;
    font-size: 12px;
  }
}
.btn {
  display: inline-block;
  margin-bottom: 0;
  font-weight: normal;
  text-align: center;
  vertical-align: middle;
  touch-action: manipulation;
  cursor: pointer;
  background-image: none;
  border: 1px solid transparent;
  white-space: nowrap;
  padding: 6px 12px;
  font-size: 13px;
  line-height: 1.42857143;
  border-radius: 2px;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}
.btn:focus,
.btn:active:focus,
.btn.active:focus,
.btn.focus,
.btn:active.focus,
.btn.active.focus {
  outline: 5px auto -webkit-focus-ring-color;
  outline-offset: -2px;
}
.btn:hover,
.btn:focus,
.btn.focus {
  color: #333;
  text-decoration: none;
}
.btn:active,
.btn.active {
  outline: 0;
  background-image: none;
  -webkit-box-shadow: inset 0 3px 5px rgba(0, 0, 0, 0.125);
  box-shadow: inset 0 3px 5px rgba(0, 0, 0, 0.125);
}
.btn.disabled,
.btn[disabled],
fieldset[disabled] .btn {
  cursor: not-allowed;
  opacity: 0.65;
  filter: alpha(opacity=65);
  -webkit-box-shadow: none;
  box-shadow: none;
}
a.btn.disabled,
fieldset[disabled] a.btn {
  pointer-events: none;
}
.btn-default {
  color: #333;
  background-color: #fff;
  border-color: #ccc;
}
.btn-default:focus,
.btn-default.focus {
  color: #333;
  background-color: #e6e6e6;
  border-color: #8c8c8c;
}
.btn-default:hover {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
.btn-default:active,
.btn-default.active,
.open > .dropdown-toggle.btn-default {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
.btn-default:active:hover,
.btn-default.active:hover,
.open > .dropdown-toggle.btn-default:hover,
.btn-default:active:focus,
.btn-default.active:focus,
.open > .dropdown-toggle.btn-default:focus,
.btn-default:active.focus,
.btn-default.active.focus,
.open > .dropdown-toggle.btn-default.focus {
  color: #333;
  background-color: #d4d4d4;
  border-color: #8c8c8c;
}
.btn-default:active,
.btn-default.active,
.open > .dropdown-toggle.btn-default {
  background-image: none;
}
.btn-default.disabled:hover,
.btn-default[disabled]:hover,
fieldset[disabled] .btn-default:hover,
.btn-default.disabled:focus,
.btn-default[disabled]:focus,
fieldset[disabled] .btn-default:focus,
.btn-default.disabled.focus,
.btn-default[disabled].focus,
fieldset[disabled] .btn-default.focus {
  background-color: #fff;
  border-color: #ccc;
}
.btn-default .badge {
  color: #fff;
  background-color: #333;
}
.btn-primary {
  color: #fff;
  background-color: #337ab7;
  border-color: #2e6da4;
}
.btn-primary:focus,
.btn-primary.focus {
  color: #fff;
  background-color: #286090;
  border-color: #122b40;
}
.btn-primary:hover {
  color: #fff;
  background-color: #286090;
  border-color: #204d74;
}
.btn-primary:active,
.btn-primary.active,
.open > .dropdown-toggle.btn-primary {
  color: #fff;
  background-color: #286090;
  border-color: #204d74;
}
.btn-primary:active:hover,
.btn-primary.active:hover,
.open > .dropdown-toggle.btn-primary:hover,
.btn-primary:active:focus,
.btn-primary.active:focus,
.open > .dropdown-toggle.btn-primary:focus,
.btn-primary:active.focus,
.btn-primary.active.focus,
.open > .dropdown-toggle.btn-primary.focus {
  color: #fff;
  background-color: #204d74;
  border-color: #122b40;
}
.btn-primary:active,
.btn-primary.active,
.open > .dropdown-toggle.btn-primary {
  background-image: none;
}
.btn-primary.disabled:hover,
.btn-primary[disabled]:hover,
fieldset[disabled] .btn-primary:hover,
.btn-primary.disabled:focus,
.btn-primary[disabled]:focus,
fieldset[disabled] .btn-primary:focus,
.btn-primary.disabled.focus,
.btn-primary[disabled].focus,
fieldset[disabled] .btn-primary.focus {
  background-color: #337ab7;
  border-color: #2e6da4;
}
.btn-primary .badge {
  color: #337ab7;
  background-color: #fff;
}
.btn-success {
  color: #fff;
  background-color: #5cb85c;
  border-color: #4cae4c;
}
.btn-success:focus,
.btn-success.focus {
  color: #fff;
  background-color: #449d44;
  border-color: #255625;
}
.btn-success:hover {
  color: #fff;
  background-color: #449d44;
  border-color: #398439;
}
.btn-success:active,
.btn-success.active,
.open > .dropdown-toggle.btn-success {
  color: #fff;
  background-color: #449d44;
  border-color: #398439;
}
.btn-success:active:hover,
.btn-success.active:hover,
.open > .dropdown-toggle.btn-success:hover,
.btn-success:active:focus,
.btn-success.active:focus,
.open > .dropdown-toggle.btn-success:focus,
.btn-success:active.focus,
.btn-success.active.focus,
.open > .dropdown-toggle.btn-success.focus {
  color: #fff;
  background-color: #398439;
  border-color: #255625;
}
.btn-success:active,
.btn-success.active,
.open > .dropdown-toggle.btn-success {
  background-image: none;
}
.btn-success.disabled:hover,
.btn-success[disabled]:hover,
fieldset[disabled] .btn-success:hover,
.btn-success.disabled:focus,
.btn-success[disabled]:focus,
fieldset[disabled] .btn-success:focus,
.btn-success.disabled.focus,
.btn-success[disabled].focus,
fieldset[disabled] .btn-success.focus {
  background-color: #5cb85c;
  border-color: #4cae4c;
}
.btn-success .badge {
  color: #5cb85c;
  background-color: #fff;
}
.btn-info {
  color: #fff;
  background-color: #5bc0de;
  border-color: #46b8da;
}
.btn-info:focus,
.btn-info.focus {
  color: #fff;
  background-color: #31b0d5;
  border-color: #1b6d85;
}
.btn-info:hover {
  color: #fff;
  background-color: #31b0d5;
  border-color: #269abc;
}
.btn-info:active,
.btn-info.active,
.open > .dropdown-toggle.btn-info {
  color: #fff;
  background-color: #31b0d5;
  border-color: #269abc;
}
.btn-info:active:hover,
.btn-info.active:hover,
.open > .dropdown-toggle.btn-info:hover,
.btn-info:active:focus,
.btn-info.active:focus,
.open > .dropdown-toggle.btn-info:focus,
.btn-info:active.focus,
.btn-info.active.focus,
.open > .dropdown-toggle.btn-info.focus {
  color: #fff;
  background-color: #269abc;
  border-color: #1b6d85;
}
.btn-info:active,
.btn-info.active,
.open > .dropdown-toggle.btn-info {
  background-image: none;
}
.btn-info.disabled:hover,
.btn-info[disabled]:hover,
fieldset[disabled] .btn-info:hover,
.btn-info.disabled:focus,
.btn-info[disabled]:focus,
fieldset[disabled] .btn-info:focus,
.btn-info.disabled.focus,
.btn-info[disabled].focus,
fieldset[disabled] .btn-info.focus {
  background-color: #5bc0de;
  border-color: #46b8da;
}
.btn-info .badge {
  color: #5bc0de;
  background-color: #fff;
}
.btn-warning {
  color: #fff;
  background-color: #f0ad4e;
  border-color: #eea236;
}
.btn-warning:focus,
.btn-warning.focus {
  color: #fff;
  background-color: #ec971f;
  border-color: #985f0d;
}
.btn-warning:hover {
  color: #fff;
  background-color: #ec971f;
  border-color: #d58512;
}
.btn-warning:active,
.btn-warning.active,
.open > .dropdown-toggle.btn-warning {
  color: #fff;
  background-color: #ec971f;
  border-color: #d58512;
}
.btn-warning:active:hover,
.btn-warning.active:hover,
.open > .dropdown-toggle.btn-warning:hover,
.btn-warning:active:focus,
.btn-warning.active:focus,
.open > .dropdown-toggle.btn-warning:focus,
.btn-warning:active.focus,
.btn-warning.active.focus,
.open > .dropdown-toggle.btn-warning.focus {
  color: #fff;
  background-color: #d58512;
  border-color: #985f0d;
}
.btn-warning:active,
.btn-warning.active,
.open > .dropdown-toggle.btn-warning {
  background-image: none;
}
.btn-warning.disabled:hover,
.btn-warning[disabled]:hover,
fieldset[disabled] .btn-warning:hover,
.btn-warning.disabled:focus,
.btn-warning[disabled]:focus,
fieldset[disabled] .btn-warning:focus,
.btn-warning.disabled.focus,
.btn-warning[disabled].focus,
fieldset[disabled] .btn-warning.focus {
  background-color: #f0ad4e;
  border-color: #eea236;
}
.btn-warning .badge {
  color: #f0ad4e;
  background-color: #fff;
}
.btn-danger {
  color: #fff;
  background-color: #d9534f;
  border-color: #d43f3a;
}
.btn-danger:focus,
.btn-danger.focus {
  color: #fff;
  background-color: #c9302c;
  border-color: #761c19;
}
.btn-danger:hover {
  color: #fff;
  background-color: #c9302c;
  border-color: #ac2925;
}
.btn-danger:active,
.btn-danger.active,
.open > .dropdown-toggle.btn-danger {
  color: #fff;
  background-color: #c9302c;
  border-color: #ac2925;
}
.btn-danger:active:hover,
.btn-danger.active:hover,
.open > .dropdown-toggle.btn-danger:hover,
.btn-danger:active:focus,
.btn-danger.active:focus,
.open > .dropdown-toggle.btn-danger:focus,
.btn-danger:active.focus,
.btn-danger.active.focus,
.open > .dropdown-toggle.btn-danger.focus {
  color: #fff;
  background-color: #ac2925;
  border-color: #761c19;
}
.btn-danger:active,
.btn-danger.active,
.open > .dropdown-toggle.btn-danger {
  background-image: none;
}
.btn-danger.disabled:hover,
.btn-danger[disabled]:hover,
fieldset[disabled] .btn-danger:hover,
.btn-danger.disabled:focus,
.btn-danger[disabled]:focus,
fieldset[disabled] .btn-danger:focus,
.btn-danger.disabled.focus,
.btn-danger[disabled].focus,
fieldset[disabled] .btn-danger.focus {
  background-color: #d9534f;
  border-color: #d43f3a;
}
.btn-danger .badge {
  color: #d9534f;
  background-color: #fff;
}
.btn-link {
  color: #337ab7;
  font-weight: normal;
  border-radius: 0;
}
.btn-link,
.btn-link:active,
.btn-link.active,
.btn-link[disabled],
fieldset[disabled] .btn-link {
  background-color: transparent;
  -webkit-box-shadow: none;
  box-shadow: none;
}
.btn-link,
.btn-link:hover,
.btn-link:focus,
.btn-link:active {
  border-color: transparent;
}
.btn-link:hover,
.btn-link:focus {
  color: #23527c;
  text-decoration: underline;
  background-color: transparent;
}
.btn-link[disabled]:hover,
fieldset[disabled] .btn-link:hover,
.btn-link[disabled]:focus,
fieldset[disabled] .btn-link:focus {
  color: #777777;
  text-decoration: none;
}
.btn-lg,
.btn-group-lg > .btn {
  padding: 10px 16px;
  font-size: 17px;
  line-height: 1.3333333;
  border-radius: 3px;
}
.btn-sm,
.btn-group-sm > .btn {
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
}
.btn-xs,
.btn-group-xs > .btn {
  padding: 1px 5px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
}
.btn-block {
  display: block;
  width: 100%;
}
.btn-block + .btn-block {
  margin-top: 5px;
}
input[type="submit"].btn-block,
input[type="reset"].btn-block,
input[type="button"].btn-block {
  width: 100%;
}
.fade {
  opacity: 0;
  -webkit-transition: opacity 0.15s linear;
  -o-transition: opacity 0.15s linear;
  transition: opacity 0.15s linear;
}
.fade.in {
  opacity: 1;
}
.collapse {
  display: none;
}
.collapse.in {
  display: block;
}
tr.collapse.in {
  display: table-row;
}
tbody.collapse.in {
  display: table-row-group;
}
.collapsing {
  position: relative;
  height: 0;
  overflow: hidden;
  -webkit-transition-property: height, visibility;
  transition-property: height, visibility;
  -webkit-transition-duration: 0.35s;
  transition-duration: 0.35s;
  -webkit-transition-timing-function: ease;
  transition-timing-function: ease;
}
.caret {
  display: inline-block;
  width: 0;
  height: 0;
  margin-left: 2px;
  vertical-align: middle;
  border-top: 4px dashed;
  border-top: 4px solid \9;
  border-right: 4px solid transparent;
  border-left: 4px solid transparent;
}
.dropup,
.dropdown {
  position: relative;
}
.dropdown-toggle:focus {
  outline: 0;
}
.dropdown-menu {
  position: absolute;
  top: 100%;
  left: 0;
  z-index: 1000;
  display: none;
  float: left;
  min-width: 160px;
  padding: 5px 0;
  margin: 2px 0 0;
  list-style: none;
  font-size: 13px;
  text-align: left;
  background-color: #fff;
  border: 1px solid #ccc;
  border: 1px solid rgba(0, 0, 0, 0.15);
  border-radius: 2px;
  -webkit-box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);
  box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);
  background-clip: padding-box;
}
.dropdown-menu.pull-right {
  right: 0;
  left: auto;
}
.dropdown-menu .divider {
  height: 1px;
  margin: 8px 0;
  overflow: hidden;
  background-color: #e5e5e5;
}
.dropdown-menu > li > a {
  display: block;
  padding: 3px 20px;
  clear: both;
  font-weight: normal;
  line-height: 1.42857143;
  color: #333333;
  white-space: nowrap;
}
.dropdown-menu > li > a:hover,
.dropdown-menu > li > a:focus {
  text-decoration: none;
  color: #262626;
  background-color: #f5f5f5;
}
.dropdown-menu > .active > a,
.dropdown-menu > .active > a:hover,
.dropdown-menu > .active > a:focus {
  color: #fff;
  text-decoration: none;
  outline: 0;
  background-color: #337ab7;
}
.dropdown-menu > .disabled > a,
.dropdown-menu > .disabled > a:hover,
.dropdown-menu > .disabled > a:focus {
  color: #777777;
}
.dropdown-menu > .disabled > a:hover,
.dropdown-menu > .disabled > a:focus {
  text-decoration: none;
  background-color: transparent;
  background-image: none;
  filter: progid:DXImageTransform.Microsoft.gradient(enabled = false);
  cursor: not-allowed;
}
.open > .dropdown-menu {
  display: block;
}
.open > a {
  outline: 0;
}
.dropdown-menu-right {
  left: auto;
  right: 0;
}
.dropdown-menu-left {
  left: 0;
  right: auto;
}
.dropdown-header {
  display: block;
  padding: 3px 20px;
  font-size: 12px;
  line-height: 1.42857143;
  color: #777777;
  white-space: nowrap;
}
.dropdown-backdrop {
  position: fixed;
  left: 0;
  right: 0;
  bottom: 0;
  top: 0;
  z-index: 990;
}
.pull-right > .dropdown-menu {
  right: 0;
  left: auto;
}
.dropup .caret,
.navbar-fixed-bottom .dropdown .caret {
  border-top: 0;
  border-bottom: 4px dashed;
  border-bottom: 4px solid \9;
  content: "";
}
.dropup .dropdown-menu,
.navbar-fixed-bottom .dropdown .dropdown-menu {
  top: auto;
  bottom: 100%;
  margin-bottom: 2px;
}
@media (min-width: 541px) {
  .navbar-right .dropdown-menu {
    left: auto;
    right: 0;
  }
  .navbar-right .dropdown-menu-left {
    left: 0;
    right: auto;
  }
}
.btn-group,
.btn-group-vertical {
  position: relative;
  display: inline-block;
  vertical-align: middle;
}
.btn-group > .btn,
.btn-group-vertical > .btn {
  position: relative;
  float: left;
}
.btn-group > .btn:hover,
.btn-group-vertical > .btn:hover,
.btn-group > .btn:focus,
.btn-group-vertical > .btn:focus,
.btn-group > .btn:active,
.btn-group-vertical > .btn:active,
.btn-group > .btn.active,
.btn-group-vertical > .btn.active {
  z-index: 2;
}
.btn-group .btn + .btn,
.btn-group .btn + .btn-group,
.btn-group .btn-group + .btn,
.btn-group .btn-group + .btn-group {
  margin-left: -1px;
}
.btn-toolbar {
  margin-left: -5px;
}
.btn-toolbar .btn,
.btn-toolbar .btn-group,
.btn-toolbar .input-group {
  float: left;
}
.btn-toolbar > .btn,
.btn-toolbar > .btn-group,
.btn-toolbar > .input-group {
  margin-left: 5px;
}
.btn-group > .btn:not(:first-child):not(:last-child):not(.dropdown-toggle) {
  border-radius: 0;
}
.btn-group > .btn:first-child {
  margin-left: 0;
}
.btn-group > .btn:first-child:not(:last-child):not(.dropdown-toggle) {
  border-bottom-right-radius: 0;
  border-top-right-radius: 0;
}
.btn-group > .btn:last-child:not(:first-child),
.btn-group > .dropdown-toggle:not(:first-child) {
  border-bottom-left-radius: 0;
  border-top-left-radius: 0;
}
.btn-group > .btn-group {
  float: left;
}
.btn-group > .btn-group:not(:first-child):not(:last-child) > .btn {
  border-radius: 0;
}
.btn-group > .btn-group:first-child:not(:last-child) > .btn:last-child,
.btn-group > .btn-group:first-child:not(:last-child) > .dropdown-toggle {
  border-bottom-right-radius: 0;
  border-top-right-radius: 0;
}
.btn-group > .btn-group:last-child:not(:first-child) > .btn:first-child {
  border-bottom-left-radius: 0;
  border-top-left-radius: 0;
}
.btn-group .dropdown-toggle:active,
.btn-group.open .dropdown-toggle {
  outline: 0;
}
.btn-group > .btn + .dropdown-toggle {
  padding-left: 8px;
  padding-right: 8px;
}
.btn-group > .btn-lg + .dropdown-toggle {
  padding-left: 12px;
  padding-right: 12px;
}
.btn-group.open .dropdown-toggle {
  -webkit-box-shadow: inset 0 3px 5px rgba(0, 0, 0, 0.125);
  box-shadow: inset 0 3px 5px rgba(0, 0, 0, 0.125);
}
.btn-group.open .dropdown-toggle.btn-link {
  -webkit-box-shadow: none;
  box-shadow: none;
}
.btn .caret {
  margin-left: 0;
}
.btn-lg .caret {
  border-width: 5px 5px 0;
  border-bottom-width: 0;
}
.dropup .btn-lg .caret {
  border-width: 0 5px 5px;
}
.btn-group-vertical > .btn,
.btn-group-vertical > .btn-group,
.btn-group-vertical > .btn-group > .btn {
  display: block;
  float: none;
  width: 100%;
  max-width: 100%;
}
.btn-group-vertical > .btn-group > .btn {
  float: none;
}
.btn-group-vertical > .btn + .btn,
.btn-group-vertical > .btn + .btn-group,
.btn-group-vertical > .btn-group + .btn,
.btn-group-vertical > .btn-group + .btn-group {
  margin-top: -1px;
  margin-left: 0;
}
.btn-group-vertical > .btn:not(:first-child):not(:last-child) {
  border-radius: 0;
}
.btn-group-vertical > .btn:first-child:not(:last-child) {
  border-top-right-radius: 2px;
  border-top-left-radius: 2px;
  border-bottom-right-radius: 0;
  border-bottom-left-radius: 0;
}
.btn-group-vertical > .btn:last-child:not(:first-child) {
  border-top-right-radius: 0;
  border-top-left-radius: 0;
  border-bottom-right-radius: 2px;
  border-bottom-left-radius: 2px;
}
.btn-group-vertical > .btn-group:not(:first-child):not(:last-child) > .btn {
  border-radius: 0;
}
.btn-group-vertical > .btn-group:first-child:not(:last-child) > .btn:last-child,
.btn-group-vertical > .btn-group:first-child:not(:last-child) > .dropdown-toggle {
  border-bottom-right-radius: 0;
  border-bottom-left-radius: 0;
}
.btn-group-vertical > .btn-group:last-child:not(:first-child) > .btn:first-child {
  border-top-right-radius: 0;
  border-top-left-radius: 0;
}
.btn-group-justified {
  display: table;
  width: 100%;
  table-layout: fixed;
  border-collapse: separate;
}
.btn-group-justified > .btn,
.btn-group-justified > .btn-group {
  float: none;
  display: table-cell;
  width: 1%;
}
.btn-group-justified > .btn-group .btn {
  width: 100%;
}
.btn-group-justified > .btn-group .dropdown-menu {
  left: auto;
}
[data-toggle="buttons"] > .btn input[type="radio"],
[data-toggle="buttons"] > .btn-group > .btn input[type="radio"],
[data-toggle="buttons"] > .btn input[type="checkbox"],
[data-toggle="buttons"] > .btn-group > .btn input[type="checkbox"] {
  position: absolute;
  clip: rect(0, 0, 0, 0);
  pointer-events: none;
}
.input-group {
  position: relative;
  display: table;
  border-collapse: separate;
}
.input-group[class*="col-"] {
  float: none;
  padding-left: 0;
  padding-right: 0;
}
.input-group .form-control {
  position: relative;
  z-index: 2;
  float: left;
  width: 100%;
  margin-bottom: 0;
}
.input-group .form-control:focus {
  z-index: 3;
}
.input-group-lg > .form-control,
.input-group-lg > .input-group-addon,
.input-group-lg > .input-group-btn > .btn {
  height: 45px;
  padding: 10px 16px;
  font-size: 17px;
  line-height: 1.3333333;
  border-radius: 3px;
}
select.input-group-lg > .form-control,
select.input-group-lg > .input-group-addon,
select.input-group-lg > .input-group-btn > .btn {
  height: 45px;
  line-height: 45px;
}
textarea.input-group-lg > .form-control,
textarea.input-group-lg > .input-group-addon,
textarea.input-group-lg > .input-group-btn > .btn,
select[multiple].input-group-lg > .form-control,
select[multiple].input-group-lg > .input-group-addon,
select[multiple].input-group-lg > .input-group-btn > .btn {
  height: auto;
}
.input-group-sm > .form-control,
.input-group-sm > .input-group-addon,
.input-group-sm > .input-group-btn > .btn {
  height: 30px;
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
}
select.input-group-sm > .form-control,
select.input-group-sm > .input-group-addon,
select.input-group-sm > .input-group-btn > .btn {
  height: 30px;
  line-height: 30px;
}
textarea.input-group-sm > .form-control,
textarea.input-group-sm > .input-group-addon,
textarea.input-group-sm > .input-group-btn > .btn,
select[multiple].input-group-sm > .form-control,
select[multiple].input-group-sm > .input-group-addon,
select[multiple].input-group-sm > .input-group-btn > .btn {
  height: auto;
}
.input-group-addon,
.input-group-btn,
.input-group .form-control {
  display: table-cell;
}
.input-group-addon:not(:first-child):not(:last-child),
.input-group-btn:not(:first-child):not(:last-child),
.input-group .form-control:not(:first-child):not(:last-child) {
  border-radius: 0;
}
.input-group-addon,
.input-group-btn {
  width: 1%;
  white-space: nowrap;
  vertical-align: middle;
}
.input-group-addon {
  padding: 6px 12px;
  font-size: 13px;
  font-weight: normal;
  line-height: 1;
  color: #555555;
  text-align: center;
  background-color: #eeeeee;
  border: 1px solid #ccc;
  border-radius: 2px;
}
.input-group-addon.input-sm {
  padding: 5px 10px;
  font-size: 12px;
  border-radius: 1px;
}
.input-group-addon.input-lg {
  padding: 10px 16px;
  font-size: 17px;
  border-radius: 3px;
}
.input-group-addon input[type="radio"],
.input-group-addon input[type="checkbox"] {
  margin-top: 0;
}
.input-group .form-control:first-child,
.input-group-addon:first-child,
.input-group-btn:first-child > .btn,
.input-group-btn:first-child > .btn-group > .btn,
.input-group-btn:first-child > .dropdown-toggle,
.input-group-btn:last-child > .btn:not(:last-child):not(.dropdown-toggle),
.input-group-btn:last-child > .btn-group:not(:last-child) > .btn {
  border-bottom-right-radius: 0;
  border-top-right-radius: 0;
}
.input-group-addon:first-child {
  border-right: 0;
}
.input-group .form-control:last-child,
.input-group-addon:last-child,
.input-group-btn:last-child > .btn,
.input-group-btn:last-child > .btn-group > .btn,
.input-group-btn:last-child > .dropdown-toggle,
.input-group-btn:first-child > .btn:not(:first-child),
.input-group-btn:first-child > .btn-group:not(:first-child) > .btn {
  border-bottom-left-radius: 0;
  border-top-left-radius: 0;
}
.input-group-addon:last-child {
  border-left: 0;
}
.input-group-btn {
  position: relative;
  font-size: 0;
  white-space: nowrap;
}
.input-group-btn > .btn {
  position: relative;
}
.input-group-btn > .btn + .btn {
  margin-left: -1px;
}
.input-group-btn > .btn:hover,
.input-group-btn > .btn:focus,
.input-group-btn > .btn:active {
  z-index: 2;
}
.input-group-btn:first-child > .btn,
.input-group-btn:first-child > .btn-group {
  margin-right: -1px;
}
.input-group-btn:last-child > .btn,
.input-group-btn:last-child > .btn-group {
  z-index: 2;
  margin-left: -1px;
}
.nav {
  margin-bottom: 0;
  padding-left: 0;
  list-style: none;
}
.nav > li {
  position: relative;
  display: block;
}
.nav > li > a {
  position: relative;
  display: block;
  padding: 10px 15px;
}
.nav > li > a:hover,
.nav > li > a:focus {
  text-decoration: none;
  background-color: #eeeeee;
}
.nav > li.disabled > a {
  color: #777777;
}
.nav > li.disabled > a:hover,
.nav > li.disabled > a:focus {
  color: #777777;
  text-decoration: none;
  background-color: transparent;
  cursor: not-allowed;
}
.nav .open > a,
.nav .open > a:hover,
.nav .open > a:focus {
  background-color: #eeeeee;
  border-color: #337ab7;
}
.nav .nav-divider {
  height: 1px;
  margin: 8px 0;
  overflow: hidden;
  background-color: #e5e5e5;
}
.nav > li > a > img {
  max-width: none;
}
.nav-tabs {
  border-bottom: 1px solid #ddd;
}
.nav-tabs > li {
  float: left;
  margin-bottom: -1px;
}
.nav-tabs > li > a {
  margin-right: 2px;
  line-height: 1.42857143;
  border: 1px solid transparent;
  border-radius: 2px 2px 0 0;
}
.nav-tabs > li > a:hover {
  border-color: #eeeeee #eeeeee #ddd;
}
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover,
.nav-tabs > li.active > a:focus {
  color: #555555;
  background-color: #fff;
  border: 1px solid #ddd;
  border-bottom-color: transparent;
  cursor: default;
}
.nav-tabs.nav-justified {
  width: 100%;
  border-bottom: 0;
}
.nav-tabs.nav-justified > li {
  float: none;
}
.nav-tabs.nav-justified > li > a {
  text-align: center;
  margin-bottom: 5px;
}
.nav-tabs.nav-justified > .dropdown .dropdown-menu {
  top: auto;
  left: auto;
}
@media (min-width: 768px) {
  .nav-tabs.nav-justified > li {
    display: table-cell;
    width: 1%;
  }
  .nav-tabs.nav-justified > li > a {
    margin-bottom: 0;
  }
}
.nav-tabs.nav-justified > li > a {
  margin-right: 0;
  border-radius: 2px;
}
.nav-tabs.nav-justified > .active > a,
.nav-tabs.nav-justified > .active > a:hover,
.nav-tabs.nav-justified > .active > a:focus {
  border: 1px solid #ddd;
}
@media (min-width: 768px) {
  .nav-tabs.nav-justified > li > a {
    border-bottom: 1px solid #ddd;
    border-radius: 2px 2px 0 0;
  }
  .nav-tabs.nav-justified > .active > a,
  .nav-tabs.nav-justified > .active > a:hover,
  .nav-tabs.nav-justified > .active > a:focus {
    border-bottom-color: #fff;
  }
}
.nav-pills > li {
  float: left;
}
.nav-pills > li > a {
  border-radius: 2px;
}
.nav-pills > li + li {
  margin-left: 2px;
}
.nav-pills > li.active > a,
.nav-pills > li.active > a:hover,
.nav-pills > li.active > a:focus {
  color: #fff;
  background-color: #337ab7;
}
.nav-stacked > li {
  float: none;
}
.nav-stacked > li + li {
  margin-top: 2px;
  margin-left: 0;
}
.nav-justified {
  width: 100%;
}
.nav-justified > li {
  float: none;
}
.nav-justified > li > a {
  text-align: center;
  margin-bottom: 5px;
}
.nav-justified > .dropdown .dropdown-menu {
  top: auto;
  left: auto;
}
@media (min-width: 768px) {
  .nav-justified > li {
    display: table-cell;
    width: 1%;
  }
  .nav-justified > li > a {
    margin-bottom: 0;
  }
}
.nav-tabs-justified {
  border-bottom: 0;
}
.nav-tabs-justified > li > a {
  margin-right: 0;
  border-radius: 2px;
}
.nav-tabs-justified > .active > a,
.nav-tabs-justified > .active > a:hover,
.nav-tabs-justified > .active > a:focus {
  border: 1px solid #ddd;
}
@media (min-width: 768px) {
  .nav-tabs-justified > li > a {
    border-bottom: 1px solid #ddd;
    border-radius: 2px 2px 0 0;
  }
  .nav-tabs-justified > .active > a,
  .nav-tabs-justified > .active > a:hover,
  .nav-tabs-justified > .active > a:focus {
    border-bottom-color: #fff;
  }
}
.tab-content > .tab-pane {
  display: none;
}
.tab-content > .active {
  display: block;
}
.nav-tabs .dropdown-menu {
  margin-top: -1px;
  border-top-right-radius: 0;
  border-top-left-radius: 0;
}
.navbar {
  position: relative;
  min-height: 30px;
  margin-bottom: 18px;
  border: 1px solid transparent;
}
@media (min-width: 541px) {
  .navbar {
    border-radius: 2px;
  }
}
@media (min-width: 541px) {
  .navbar-header {
    float: left;
  }
}
.navbar-collapse {
  overflow-x: visible;
  padding-right: 0px;
  padding-left: 0px;
  border-top: 1px solid transparent;
  box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.1);
  -webkit-overflow-scrolling: touch;
}
.navbar-collapse.in {
  overflow-y: auto;
}
@media (min-width: 541px) {
  .navbar-collapse {
    width: auto;
    border-top: 0;
    box-shadow: none;
  }
  .navbar-collapse.collapse {
    display: block !important;
    height: auto !important;
    padding-bottom: 0;
    overflow: visible !important;
  }
  .navbar-collapse.in {
    overflow-y: visible;
  }
  .navbar-fixed-top .navbar-collapse,
  .navbar-static-top .navbar-collapse,
  .navbar-fixed-bottom .navbar-collapse {
    padding-left: 0;
    padding-right: 0;
  }
}
.navbar-fixed-top .navbar-collapse,
.navbar-fixed-bottom .navbar-collapse {
  max-height: 340px;
}
@media (max-device-width: 540px) and (orientation: landscape) {
  .navbar-fixed-top .navbar-collapse,
  .navbar-fixed-bottom .navbar-collapse {
    max-height: 200px;
  }
}
.container > .navbar-header,
.container-fluid > .navbar-header,
.container > .navbar-collapse,
.container-fluid > .navbar-collapse {
  margin-right: 0px;
  margin-left: 0px;
}
@media (min-width: 541px) {
  .container > .navbar-header,
  .container-fluid > .navbar-header,
  .container > .navbar-collapse,
  .container-fluid > .navbar-collapse {
    margin-right: 0;
    margin-left: 0;
  }
}
.navbar-static-top {
  z-index: 1000;
  border-width: 0 0 1px;
}
@media (min-width: 541px) {
  .navbar-static-top {
    border-radius: 0;
  }
}
.navbar-fixed-top,
.navbar-fixed-bottom {
  position: fixed;
  right: 0;
  left: 0;
  z-index: 1030;
}
@media (min-width: 541px) {
  .navbar-fixed-top,
  .navbar-fixed-bottom {
    border-radius: 0;
  }
}
.navbar-fixed-top {
  top: 0;
  border-width: 0 0 1px;
}
.navbar-fixed-bottom {
  bottom: 0;
  margin-bottom: 0;
  border-width: 1px 0 0;
}
.navbar-brand {
  float: left;
  padding: 6px 0px;
  font-size: 17px;
  line-height: 18px;
  height: 30px;
}
.navbar-brand:hover,
.navbar-brand:focus {
  text-decoration: none;
}
.navbar-brand > img {
  display: block;
}
@media (min-width: 541px) {
  .navbar > .container .navbar-brand,
  .navbar > .container-fluid .navbar-brand {
    margin-left: 0px;
  }
}
.navbar-toggle {
  position: relative;
  float: right;
  margin-right: 0px;
  padding: 9px 10px;
  margin-top: -2px;
  margin-bottom: -2px;
  background-color: transparent;
  background-image: none;
  border: 1px solid transparent;
  border-radius: 2px;
}
.navbar-toggle:focus {
  outline: 0;
}
.navbar-toggle .icon-bar {
  display: block;
  width: 22px;
  height: 2px;
  border-radius: 1px;
}
.navbar-toggle .icon-bar + .icon-bar {
  margin-top: 4px;
}
@media (min-width: 541px) {
  .navbar-toggle {
    display: none;
  }
}
.navbar-nav {
  margin: 3px 0px;
}
.navbar-nav > li > a {
  padding-top: 10px;
  padding-bottom: 10px;
  line-height: 18px;
}
@media (max-width: 540px) {
  .navbar-nav .open .dropdown-menu {
    position: static;
    float: none;
    width: auto;
    margin-top: 0;
    background-color: transparent;
    border: 0;
    box-shadow: none;
  }
  .navbar-nav .open .dropdown-menu > li > a,
  .navbar-nav .open .dropdown-menu .dropdown-header {
    padding: 5px 15px 5px 25px;
  }
  .navbar-nav .open .dropdown-menu > li > a {
    line-height: 18px;
  }
  .navbar-nav .open .dropdown-menu > li > a:hover,
  .navbar-nav .open .dropdown-menu > li > a:focus {
    background-image: none;
  }
}
@media (min-width: 541px) {
  .navbar-nav {
    float: left;
    margin: 0;
  }
  .navbar-nav > li {
    float: left;
  }
  .navbar-nav > li > a {
    padding-top: 6px;
    padding-bottom: 6px;
  }
}
.navbar-form {
  margin-left: 0px;
  margin-right: 0px;
  padding: 10px 0px;
  border-top: 1px solid transparent;
  border-bottom: 1px solid transparent;
  -webkit-box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.1), 0 1px 0 rgba(255, 255, 255, 0.1);
  box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.1), 0 1px 0 rgba(255, 255, 255, 0.1);
  margin-top: -1px;
  margin-bottom: -1px;
}
@media (min-width: 768px) {
  .navbar-form .form-group {
    display: inline-block;
    margin-bottom: 0;
    vertical-align: middle;
  }
  .navbar-form .form-control {
    display: inline-block;
    width: auto;
    vertical-align: middle;
  }
  .navbar-form .form-control-static {
    display: inline-block;
  }
  .navbar-form .input-group {
    display: inline-table;
    vertical-align: middle;
  }
  .navbar-form .input-group .input-group-addon,
  .navbar-form .input-group .input-group-btn,
  .navbar-form .input-group .form-control {
    width: auto;
  }
  .navbar-form .input-group > .form-control {
    width: 100%;
  }
  .navbar-form .control-label {
    margin-bottom: 0;
    vertical-align: middle;
  }
  .navbar-form .radio,
  .navbar-form .checkbox {
    display: inline-block;
    margin-top: 0;
    margin-bottom: 0;
    vertical-align: middle;
  }
  .navbar-form .radio label,
  .navbar-form .checkbox label {
    padding-left: 0;
  }
  .navbar-form .radio input[type="radio"],
  .navbar-form .checkbox input[type="checkbox"] {
    position: relative;
    margin-left: 0;
  }
  .navbar-form .has-feedback .form-control-feedback {
    top: 0;
  }
}
@media (max-width: 540px) {
  .navbar-form .form-group {
    margin-bottom: 5px;
  }
  .navbar-form .form-group:last-child {
    margin-bottom: 0;
  }
}
@media (min-width: 541px) {
  .navbar-form {
    width: auto;
    border: 0;
    margin-left: 0;
    margin-right: 0;
    padding-top: 0;
    padding-bottom: 0;
    -webkit-box-shadow: none;
    box-shadow: none;
  }
}
.navbar-nav > li > .dropdown-menu {
  margin-top: 0;
  border-top-right-radius: 0;
  border-top-left-radius: 0;
}
.navbar-fixed-bottom .navbar-nav > li > .dropdown-menu {
  margin-bottom: 0;
  border-top-right-radius: 2px;
  border-top-left-radius: 2px;
  border-bottom-right-radius: 0;
  border-bottom-left-radius: 0;
}
.navbar-btn {
  margin-top: -1px;
  margin-bottom: -1px;
}
.navbar-btn.btn-sm {
  margin-top: 0px;
  margin-bottom: 0px;
}
.navbar-btn.btn-xs {
  margin-top: 4px;
  margin-bottom: 4px;
}
.navbar-text {
  margin-top: 6px;
  margin-bottom: 6px;
}
@media (min-width: 541px) {
  .navbar-text {
    float: left;
    margin-left: 0px;
    margin-right: 0px;
  }
}
@media (min-width: 541px) {
  .navbar-left {
    float: left !important;
    float: left;
  }
  .navbar-right {
    float: right !important;
    float: right;
    margin-right: 0px;
  }
  .navbar-right ~ .navbar-right {
    margin-right: 0;
  }
}
.navbar-default {
  background-color: #f8f8f8;
  border-color: #e7e7e7;
}
.navbar-default .navbar-brand {
  color: #777;
}
.navbar-default .navbar-brand:hover,
.navbar-default .navbar-brand:focus {
  color: #5e5e5e;
  background-color: transparent;
}
.navbar-default .navbar-text {
  color: #777;
}
.navbar-default .navbar-nav > li > a {
  color: #777;
}
.navbar-default .navbar-nav > li > a:hover,
.navbar-default .navbar-nav > li > a:focus {
  color: #333;
  background-color: transparent;
}
.navbar-default .navbar-nav > .active > a,
.navbar-default .navbar-nav > .active > a:hover,
.navbar-default .navbar-nav > .active > a:focus {
  color: #555;
  background-color: #e7e7e7;
}
.navbar-default .navbar-nav > .disabled > a,
.navbar-default .navbar-nav > .disabled > a:hover,
.navbar-default .navbar-nav > .disabled > a:focus {
  color: #ccc;
  background-color: transparent;
}
.navbar-default .navbar-toggle {
  border-color: #ddd;
}
.navbar-default .navbar-toggle:hover,
.navbar-default .navbar-toggle:focus {
  background-color: #ddd;
}
.navbar-default .navbar-toggle .icon-bar {
  background-color: #888;
}
.navbar-default .navbar-collapse,
.navbar-default .navbar-form {
  border-color: #e7e7e7;
}
.navbar-default .navbar-nav > .open > a,
.navbar-default .navbar-nav > .open > a:hover,
.navbar-default .navbar-nav > .open > a:focus {
  background-color: #e7e7e7;
  color: #555;
}
@media (max-width: 540px) {
  .navbar-default .navbar-nav .open .dropdown-menu > li > a {
    color: #777;
  }
  .navbar-default .navbar-nav .open .dropdown-menu > li > a:hover,
  .navbar-default .navbar-nav .open .dropdown-menu > li > a:focus {
    color: #333;
    background-color: transparent;
  }
  .navbar-default .navbar-nav .open .dropdown-menu > .active > a,
  .navbar-default .navbar-nav .open .dropdown-menu > .active > a:hover,
  .navbar-default .navbar-nav .open .dropdown-menu > .active > a:focus {
    color: #555;
    background-color: #e7e7e7;
  }
  .navbar-default .navbar-nav .open .dropdown-menu > .disabled > a,
  .navbar-default .navbar-nav .open .dropdown-menu > .disabled > a:hover,
  .navbar-default .navbar-nav .open .dropdown-menu > .disabled > a:focus {
    color: #ccc;
    background-color: transparent;
  }
}
.navbar-default .navbar-link {
  color: #777;
}
.navbar-default .navbar-link:hover {
  color: #333;
}
.navbar-default .btn-link {
  color: #777;
}
.navbar-default .btn-link:hover,
.navbar-default .btn-link:focus {
  color: #333;
}
.navbar-default .btn-link[disabled]:hover,
fieldset[disabled] .navbar-default .btn-link:hover,
.navbar-default .btn-link[disabled]:focus,
fieldset[disabled] .navbar-default .btn-link:focus {
  color: #ccc;
}
.navbar-inverse {
  background-color: #222;
  border-color: #080808;
}
.navbar-inverse .navbar-brand {
  color: #9d9d9d;
}
.navbar-inverse .navbar-brand:hover,
.navbar-inverse .navbar-brand:focus {
  color: #fff;
  background-color: transparent;
}
.navbar-inverse .navbar-text {
  color: #9d9d9d;
}
.navbar-inverse .navbar-nav > li > a {
  color: #9d9d9d;
}
.navbar-inverse .navbar-nav > li > a:hover,
.navbar-inverse .navbar-nav > li > a:focus {
  color: #fff;
  background-color: transparent;
}
.navbar-inverse .navbar-nav > .active > a,
.navbar-inverse .navbar-nav > .active > a:hover,
.navbar-inverse .navbar-nav > .active > a:focus {
  color: #fff;
  background-color: #080808;
}
.navbar-inverse .navbar-nav > .disabled > a,
.navbar-inverse .navbar-nav > .disabled > a:hover,
.navbar-inverse .navbar-nav > .disabled > a:focus {
  color: #444;
  background-color: transparent;
}
.navbar-inverse .navbar-toggle {
  border-color: #333;
}
.navbar-inverse .navbar-toggle:hover,
.navbar-inverse .navbar-toggle:focus {
  background-color: #333;
}
.navbar-inverse .navbar-toggle .icon-bar {
  background-color: #fff;
}
.navbar-inverse .navbar-collapse,
.navbar-inverse .navbar-form {
  border-color: #101010;
}
.navbar-inverse .navbar-nav > .open > a,
.navbar-inverse .navbar-nav > .open > a:hover,
.navbar-inverse .navbar-nav > .open > a:focus {
  background-color: #080808;
  color: #fff;
}
@media (max-width: 540px) {
  .navbar-inverse .navbar-nav .open .dropdown-menu > .dropdown-header {
    border-color: #080808;
  }
  .navbar-inverse .navbar-nav .open .dropdown-menu .divider {
    background-color: #080808;
  }
  .navbar-inverse .navbar-nav .open .dropdown-menu > li > a {
    color: #9d9d9d;
  }
  .navbar-inverse .navbar-nav .open .dropdown-menu > li > a:hover,
  .navbar-inverse .navbar-nav .open .dropdown-menu > li > a:focus {
    color: #fff;
    background-color: transparent;
  }
  .navbar-inverse .navbar-nav .open .dropdown-menu > .active > a,
  .navbar-inverse .navbar-nav .open .dropdown-menu > .active > a:hover,
  .navbar-inverse .navbar-nav .open .dropdown-menu > .active > a:focus {
    color: #fff;
    background-color: #080808;
  }
  .navbar-inverse .navbar-nav .open .dropdown-menu > .disabled > a,
  .navbar-inverse .navbar-nav .open .dropdown-menu > .disabled > a:hover,
  .navbar-inverse .navbar-nav .open .dropdown-menu > .disabled > a:focus {
    color: #444;
    background-color: transparent;
  }
}
.navbar-inverse .navbar-link {
  color: #9d9d9d;
}
.navbar-inverse .navbar-link:hover {
  color: #fff;
}
.navbar-inverse .btn-link {
  color: #9d9d9d;
}
.navbar-inverse .btn-link:hover,
.navbar-inverse .btn-link:focus {
  color: #fff;
}
.navbar-inverse .btn-link[disabled]:hover,
fieldset[disabled] .navbar-inverse .btn-link:hover,
.navbar-inverse .btn-link[disabled]:focus,
fieldset[disabled] .navbar-inverse .btn-link:focus {
  color: #444;
}
.breadcrumb {
  padding: 8px 15px;
  margin-bottom: 18px;
  list-style: none;
  background-color: #f5f5f5;
  border-radius: 2px;
}
.breadcrumb > li {
  display: inline-block;
}
.breadcrumb > li + li:before {
  content: "/\00a0";
  padding: 0 5px;
  color: #5e5e5e;
}
.breadcrumb > .active {
  color: #777777;
}
.pagination {
  display: inline-block;
  padding-left: 0;
  margin: 18px 0;
  border-radius: 2px;
}
.pagination > li {
  display: inline;
}
.pagination > li > a,
.pagination > li > span {
  position: relative;
  float: left;
  padding: 6px 12px;
  line-height: 1.42857143;
  text-decoration: none;
  color: #337ab7;
  background-color: #fff;
  border: 1px solid #ddd;
  margin-left: -1px;
}
.pagination > li:first-child > a,
.pagination > li:first-child > span {
  margin-left: 0;
  border-bottom-left-radius: 2px;
  border-top-left-radius: 2px;
}
.pagination > li:last-child > a,
.pagination > li:last-child > span {
  border-bottom-right-radius: 2px;
  border-top-right-radius: 2px;
}
.pagination > li > a:hover,
.pagination > li > span:hover,
.pagination > li > a:focus,
.pagination > li > span:focus {
  z-index: 2;
  color: #23527c;
  background-color: #eeeeee;
  border-color: #ddd;
}
.pagination > .active > a,
.pagination > .active > span,
.pagination > .active > a:hover,
.pagination > .active > span:hover,
.pagination > .active > a:focus,
.pagination > .active > span:focus {
  z-index: 3;
  color: #fff;
  background-color: #337ab7;
  border-color: #337ab7;
  cursor: default;
}
.pagination > .disabled > span,
.pagination > .disabled > span:hover,
.pagination > .disabled > span:focus,
.pagination > .disabled > a,
.pagination > .disabled > a:hover,
.pagination > .disabled > a:focus {
  color: #777777;
  background-color: #fff;
  border-color: #ddd;
  cursor: not-allowed;
}
.pagination-lg > li > a,
.pagination-lg > li > span {
  padding: 10px 16px;
  font-size: 17px;
  line-height: 1.3333333;
}
.pagination-lg > li:first-child > a,
.pagination-lg > li:first-child > span {
  border-bottom-left-radius: 3px;
  border-top-left-radius: 3px;
}
.pagination-lg > li:last-child > a,
.pagination-lg > li:last-child > span {
  border-bottom-right-radius: 3px;
  border-top-right-radius: 3px;
}
.pagination-sm > li > a,
.pagination-sm > li > span {
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
}
.pagination-sm > li:first-child > a,
.pagination-sm > li:first-child > span {
  border-bottom-left-radius: 1px;
  border-top-left-radius: 1px;
}
.pagination-sm > li:last-child > a,
.pagination-sm > li:last-child > span {
  border-bottom-right-radius: 1px;
  border-top-right-radius: 1px;
}
.pager {
  padding-left: 0;
  margin: 18px 0;
  list-style: none;
  text-align: center;
}
.pager li {
  display: inline;
}
.pager li > a,
.pager li > span {
  display: inline-block;
  padding: 5px 14px;
  background-color: #fff;
  border: 1px solid #ddd;
  border-radius: 15px;
}
.pager li > a:hover,
.pager li > a:focus {
  text-decoration: none;
  background-color: #eeeeee;
}
.pager .next > a,
.pager .next > span {
  float: right;
}
.pager .previous > a,
.pager .previous > span {
  float: left;
}
.pager .disabled > a,
.pager .disabled > a:hover,
.pager .disabled > a:focus,
.pager .disabled > span {
  color: #777777;
  background-color: #fff;
  cursor: not-allowed;
}
.label {
  display: inline;
  padding: .2em .6em .3em;
  font-size: 75%;
  font-weight: bold;
  line-height: 1;
  color: #fff;
  text-align: center;
  white-space: nowrap;
  vertical-align: baseline;
  border-radius: .25em;
}
a.label:hover,
a.label:focus {
  color: #fff;
  text-decoration: none;
  cursor: pointer;
}
.label:empty {
  display: none;
}
.btn .label {
  position: relative;
  top: -1px;
}
.label-default {
  background-color: #777777;
}
.label-default[href]:hover,
.label-default[href]:focus {
  background-color: #5e5e5e;
}
.label-primary {
  background-color: #337ab7;
}
.label-primary[href]:hover,
.label-primary[href]:focus {
  background-color: #286090;
}
.label-success {
  background-color: #5cb85c;
}
.label-success[href]:hover,
.label-success[href]:focus {
  background-color: #449d44;
}
.label-info {
  background-color: #5bc0de;
}
.label-info[href]:hover,
.label-info[href]:focus {
  background-color: #31b0d5;
}
.label-warning {
  background-color: #f0ad4e;
}
.label-warning[href]:hover,
.label-warning[href]:focus {
  background-color: #ec971f;
}
.label-danger {
  background-color: #d9534f;
}
.label-danger[href]:hover,
.label-danger[href]:focus {
  background-color: #c9302c;
}
.badge {
  display: inline-block;
  min-width: 10px;
  padding: 3px 7px;
  font-size: 12px;
  font-weight: bold;
  color: #fff;
  line-height: 1;
  vertical-align: middle;
  white-space: nowrap;
  text-align: center;
  background-color: #777777;
  border-radius: 10px;
}
.badge:empty {
  display: none;
}
.btn .badge {
  position: relative;
  top: -1px;
}
.btn-xs .badge,
.btn-group-xs > .btn .badge {
  top: 0;
  padding: 1px 5px;
}
a.badge:hover,
a.badge:focus {
  color: #fff;
  text-decoration: none;
  cursor: pointer;
}
.list-group-item.active > .badge,
.nav-pills > .active > a > .badge {
  color: #337ab7;
  background-color: #fff;
}
.list-group-item > .badge {
  float: right;
}
.list-group-item > .badge + .badge {
  margin-right: 5px;
}
.nav-pills > li > a > .badge {
  margin-left: 3px;
}
.jumbotron {
  padding-top: 30px;
  padding-bottom: 30px;
  margin-bottom: 30px;
  color: inherit;
  background-color: #eeeeee;
}
.jumbotron h1,
.jumbotron .h1 {
  color: inherit;
}
.jumbotron p {
  margin-bottom: 15px;
  font-size: 20px;
  font-weight: 200;
}
.jumbotron > hr {
  border-top-color: #d5d5d5;
}
.container .jumbotron,
.container-fluid .jumbotron {
  border-radius: 3px;
  padding-left: 0px;
  padding-right: 0px;
}
.jumbotron .container {
  max-width: 100%;
}
@media screen and (min-width: 768px) {
  .jumbotron {
    padding-top: 48px;
    padding-bottom: 48px;
  }
  .container .jumbotron,
  .container-fluid .jumbotron {
    padding-left: 60px;
    padding-right: 60px;
  }
  .jumbotron h1,
  .jumbotron .h1 {
    font-size: 59px;
  }
}
.thumbnail {
  display: block;
  padding: 4px;
  margin-bottom: 18px;
  line-height: 1.42857143;
  background-color: #fff;
  border: 1px solid #ddd;
  border-radius: 2px;
  -webkit-transition: border 0.2s ease-in-out;
  -o-transition: border 0.2s ease-in-out;
  transition: border 0.2s ease-in-out;
}
.thumbnail > img,
.thumbnail a > img {
  margin-left: auto;
  margin-right: auto;
}
a.thumbnail:hover,
a.thumbnail:focus,
a.thumbnail.active {
  border-color: #337ab7;
}
.thumbnail .caption {
  padding: 9px;
  color: #000;
}
.alert {
  padding: 15px;
  margin-bottom: 18px;
  border: 1px solid transparent;
  border-radius: 2px;
}
.alert h4 {
  margin-top: 0;
  color: inherit;
}
.alert .alert-link {
  font-weight: bold;
}
.alert > p,
.alert > ul {
  margin-bottom: 0;
}
.alert > p + p {
  margin-top: 5px;
}
.alert-dismissable,
.alert-dismissible {
  padding-right: 35px;
}
.alert-dismissable .close,
.alert-dismissible .close {
  position: relative;
  top: -2px;
  right: -21px;
  color: inherit;
}
.alert-success {
  background-color: #dff0d8;
  border-color: #d6e9c6;
  color: #3c763d;
}
.alert-success hr {
  border-top-color: #c9e2b3;
}
.alert-success .alert-link {
  color: #2b542c;
}
.alert-info {
  background-color: #d9edf7;
  border-color: #bce8f1;
  color: #31708f;
}
.alert-info hr {
  border-top-color: #a6e1ec;
}
.alert-info .alert-link {
  color: #245269;
}
.alert-warning {
  background-color: #fcf8e3;
  border-color: #faebcc;
  color: #8a6d3b;
}
.alert-warning hr {
  border-top-color: #f7e1b5;
}
.alert-warning .alert-link {
  color: #66512c;
}
.alert-danger {
  background-color: #f2dede;
  border-color: #ebccd1;
  color: #a94442;
}
.alert-danger hr {
  border-top-color: #e4b9c0;
}
.alert-danger .alert-link {
  color: #843534;
}
@-webkit-keyframes progress-bar-stripes {
  from {
    background-position: 40px 0;
  }
  to {
    background-position: 0 0;
  }
}
@keyframes progress-bar-stripes {
  from {
    background-position: 40px 0;
  }
  to {
    background-position: 0 0;
  }
}
.progress {
  overflow: hidden;
  height: 18px;
  margin-bottom: 18px;
  background-color: #f5f5f5;
  border-radius: 2px;
  -webkit-box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.1);
  box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.1);
}
.progress-bar {
  float: left;
  width: 0%;
  height: 100%;
  font-size: 12px;
  line-height: 18px;
  color: #fff;
  text-align: center;
  background-color: #337ab7;
  -webkit-box-shadow: inset 0 -1px 0 rgba(0, 0, 0, 0.15);
  box-shadow: inset 0 -1px 0 rgba(0, 0, 0, 0.15);
  -webkit-transition: width 0.6s ease;
  -o-transition: width 0.6s ease;
  transition: width 0.6s ease;
}
.progress-striped .progress-bar,
.progress-bar-striped {
  background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-size: 40px 40px;
}
.progress.active .progress-bar,
.progress-bar.active {
  -webkit-animation: progress-bar-stripes 2s linear infinite;
  -o-animation: progress-bar-stripes 2s linear infinite;
  animation: progress-bar-stripes 2s linear infinite;
}
.progress-bar-success {
  background-color: #5cb85c;
}
.progress-striped .progress-bar-success {
  background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
}
.progress-bar-info {
  background-color: #5bc0de;
}
.progress-striped .progress-bar-info {
  background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
}
.progress-bar-warning {
  background-color: #f0ad4e;
}
.progress-striped .progress-bar-warning {
  background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
}
.progress-bar-danger {
  background-color: #d9534f;
}
.progress-striped .progress-bar-danger {
  background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
}
.media {
  margin-top: 15px;
}
.media:first-child {
  margin-top: 0;
}
.media,
.media-body {
  zoom: 1;
  overflow: hidden;
}
.media-body {
  width: 10000px;
}
.media-object {
  display: block;
}
.media-object.img-thumbnail {
  max-width: none;
}
.media-right,
.media > .pull-right {
  padding-left: 10px;
}
.media-left,
.media > .pull-left {
  padding-right: 10px;
}
.media-left,
.media-right,
.media-body {
  display: table-cell;
  vertical-align: top;
}
.media-middle {
  vertical-align: middle;
}
.media-bottom {
  vertical-align: bottom;
}
.media-heading {
  margin-top: 0;
  margin-bottom: 5px;
}
.media-list {
  padding-left: 0;
  list-style: none;
}
.list-group {
  margin-bottom: 20px;
  padding-left: 0;
}
.list-group-item {
  position: relative;
  display: block;
  padding: 10px 15px;
  margin-bottom: -1px;
  background-color: #fff;
  border: 1px solid #ddd;
}
.list-group-item:first-child {
  border-top-right-radius: 2px;
  border-top-left-radius: 2px;
}
.list-group-item:last-child {
  margin-bottom: 0;
  border-bottom-right-radius: 2px;
  border-bottom-left-radius: 2px;
}
a.list-group-item,
button.list-group-item {
  color: #555;
}
a.list-group-item .list-group-item-heading,
button.list-group-item .list-group-item-heading {
  color: #333;
}
a.list-group-item:hover,
button.list-group-item:hover,
a.list-group-item:focus,
button.list-group-item:focus {
  text-decoration: none;
  color: #555;
  background-color: #f5f5f5;
}
button.list-group-item {
  width: 100%;
  text-align: left;
}
.list-group-item.disabled,
.list-group-item.disabled:hover,
.list-group-item.disabled:focus {
  background-color: #eeeeee;
  color: #777777;
  cursor: not-allowed;
}
.list-group-item.disabled .list-group-item-heading,
.list-group-item.disabled:hover .list-group-item-heading,
.list-group-item.disabled:focus .list-group-item-heading {
  color: inherit;
}
.list-group-item.disabled .list-group-item-text,
.list-group-item.disabled:hover .list-group-item-text,
.list-group-item.disabled:focus .list-group-item-text {
  color: #777777;
}
.list-group-item.active,
.list-group-item.active:hover,
.list-group-item.active:focus {
  z-index: 2;
  color: #fff;
  background-color: #337ab7;
  border-color: #337ab7;
}
.list-group-item.active .list-group-item-heading,
.list-group-item.active:hover .list-group-item-heading,
.list-group-item.active:focus .list-group-item-heading,
.list-group-item.active .list-group-item-heading > small,
.list-group-item.active:hover .list-group-item-heading > small,
.list-group-item.active:focus .list-group-item-heading > small,
.list-group-item.active .list-group-item-heading > .small,
.list-group-item.active:hover .list-group-item-heading > .small,
.list-group-item.active:focus .list-group-item-heading > .small {
  color: inherit;
}
.list-group-item.active .list-group-item-text,
.list-group-item.active:hover .list-group-item-text,
.list-group-item.active:focus .list-group-item-text {
  color: #c7ddef;
}
.list-group-item-success {
  color: #3c763d;
  background-color: #dff0d8;
}
a.list-group-item-success,
button.list-group-item-success {
  color: #3c763d;
}
a.list-group-item-success .list-group-item-heading,
button.list-group-item-success .list-group-item-heading {
  color: inherit;
}
a.list-group-item-success:hover,
button.list-group-item-success:hover,
a.list-group-item-success:focus,
button.list-group-item-success:focus {
  color: #3c763d;
  background-color: #d0e9c6;
}
a.list-group-item-success.active,
button.list-group-item-success.active,
a.list-group-item-success.active:hover,
button.list-group-item-success.active:hover,
a.list-group-item-success.active:focus,
button.list-group-item-success.active:focus {
  color: #fff;
  background-color: #3c763d;
  border-color: #3c763d;
}
.list-group-item-info {
  color: #31708f;
  background-color: #d9edf7;
}
a.list-group-item-info,
button.list-group-item-info {
  color: #31708f;
}
a.list-group-item-info .list-group-item-heading,
button.list-group-item-info .list-group-item-heading {
  color: inherit;
}
a.list-group-item-info:hover,
button.list-group-item-info:hover,
a.list-group-item-info:focus,
button.list-group-item-info:focus {
  color: #31708f;
  background-color: #c4e3f3;
}
a.list-group-item-info.active,
button.list-group-item-info.active,
a.list-group-item-info.active:hover,
button.list-group-item-info.active:hover,
a.list-group-item-info.active:focus,
button.list-group-item-info.active:focus {
  color: #fff;
  background-color: #31708f;
  border-color: #31708f;
}
.list-group-item-warning {
  color: #8a6d3b;
  background-color: #fcf8e3;
}
a.list-group-item-warning,
button.list-group-item-warning {
  color: #8a6d3b;
}
a.list-group-item-warning .list-group-item-heading,
button.list-group-item-warning .list-group-item-heading {
  color: inherit;
}
a.list-group-item-warning:hover,
button.list-group-item-warning:hover,
a.list-group-item-warning:focus,
button.list-group-item-warning:focus {
  color: #8a6d3b;
  background-color: #faf2cc;
}
a.list-group-item-warning.active,
button.list-group-item-warning.active,
a.list-group-item-warning.active:hover,
button.list-group-item-warning.active:hover,
a.list-group-item-warning.active:focus,
button.list-group-item-warning.active:focus {
  color: #fff;
  background-color: #8a6d3b;
  border-color: #8a6d3b;
}
.list-group-item-danger {
  color: #a94442;
  background-color: #f2dede;
}
a.list-group-item-danger,
button.list-group-item-danger {
  color: #a94442;
}
a.list-group-item-danger .list-group-item-heading,
button.list-group-item-danger .list-group-item-heading {
  color: inherit;
}
a.list-group-item-danger:hover,
button.list-group-item-danger:hover,
a.list-group-item-danger:focus,
button.list-group-item-danger:focus {
  color: #a94442;
  background-color: #ebcccc;
}
a.list-group-item-danger.active,
button.list-group-item-danger.active,
a.list-group-item-danger.active:hover,
button.list-group-item-danger.active:hover,
a.list-group-item-danger.active:focus,
button.list-group-item-danger.active:focus {
  color: #fff;
  background-color: #a94442;
  border-color: #a94442;
}
.list-group-item-heading {
  margin-top: 0;
  margin-bottom: 5px;
}
.list-group-item-text {
  margin-bottom: 0;
  line-height: 1.3;
}
.panel {
  margin-bottom: 18px;
  background-color: #fff;
  border: 1px solid transparent;
  border-radius: 2px;
  -webkit-box-shadow: 0 1px 1px rgba(0, 0, 0, 0.05);
  box-shadow: 0 1px 1px rgba(0, 0, 0, 0.05);
}
.panel-body {
  padding: 15px;
}
.panel-heading {
  padding: 10px 15px;
  border-bottom: 1px solid transparent;
  border-top-right-radius: 1px;
  border-top-left-radius: 1px;
}
.panel-heading > .dropdown .dropdown-toggle {
  color: inherit;
}
.panel-title {
  margin-top: 0;
  margin-bottom: 0;
  font-size: 15px;
  color: inherit;
}
.panel-title > a,
.panel-title > small,
.panel-title > .small,
.panel-title > small > a,
.panel-title > .small > a {
  color: inherit;
}
.panel-footer {
  padding: 10px 15px;
  background-color: #f5f5f5;
  border-top: 1px solid #ddd;
  border-bottom-right-radius: 1px;
  border-bottom-left-radius: 1px;
}
.panel > .list-group,
.panel > .panel-collapse > .list-group {
  margin-bottom: 0;
}
.panel > .list-group .list-group-item,
.panel > .panel-collapse > .list-group .list-group-item {
  border-width: 1px 0;
  border-radius: 0;
}
.panel > .list-group:first-child .list-group-item:first-child,
.panel > .panel-collapse > .list-group:first-child .list-group-item:first-child {
  border-top: 0;
  border-top-right-radius: 1px;
  border-top-left-radius: 1px;
}
.panel > .list-group:last-child .list-group-item:last-child,
.panel > .panel-collapse > .list-group:last-child .list-group-item:last-child {
  border-bottom: 0;
  border-bottom-right-radius: 1px;
  border-bottom-left-radius: 1px;
}
.panel > .panel-heading + .panel-collapse > .list-group .list-group-item:first-child {
  border-top-right-radius: 0;
  border-top-left-radius: 0;
}
.panel-heading + .list-group .list-group-item:first-child {
  border-top-width: 0;
}
.list-group + .panel-footer {
  border-top-width: 0;
}
.panel > .table,
.panel > .table-responsive > .table,
.panel > .panel-collapse > .table {
  margin-bottom: 0;
}
.panel > .table caption,
.panel > .table-responsive > .table caption,
.panel > .panel-collapse > .table caption {
  padding-left: 15px;
  padding-right: 15px;
}
.panel > .table:first-child,
.panel > .table-responsive:first-child > .table:first-child {
  border-top-right-radius: 1px;
  border-top-left-radius: 1px;
}
.panel > .table:first-child > thead:first-child > tr:first-child,
.panel > .table-responsive:first-child > .table:first-child > thead:first-child > tr:first-child,
.panel > .table:first-child > tbody:first-child > tr:first-child,
.panel > .table-responsive:first-child > .table:first-child > tbody:first-child > tr:first-child {
  border-top-left-radius: 1px;
  border-top-right-radius: 1px;
}
.panel > .table:first-child > thead:first-child > tr:first-child td:first-child,
.panel > .table-responsive:first-child > .table:first-child > thead:first-child > tr:first-child td:first-child,
.panel > .table:first-child > tbody:first-child > tr:first-child td:first-child,
.panel > .table-responsive:first-child > .table:first-child > tbody:first-child > tr:first-child td:first-child,
.panel > .table:first-child > thead:first-child > tr:first-child th:first-child,
.panel > .table-responsive:first-child > .table:first-child > thead:first-child > tr:first-child th:first-child,
.panel > .table:first-child > tbody:first-child > tr:first-child th:first-child,
.panel > .table-responsive:first-child > .table:first-child > tbody:first-child > tr:first-child th:first-child {
  border-top-left-radius: 1px;
}
.panel > .table:first-child > thead:first-child > tr:first-child td:last-child,
.panel > .table-responsive:first-child > .table:first-child > thead:first-child > tr:first-child td:last-child,
.panel > .table:first-child > tbody:first-child > tr:first-child td:last-child,
.panel > .table-responsive:first-child > .table:first-child > tbody:first-child > tr:first-child td:last-child,
.panel > .table:first-child > thead:first-child > tr:first-child th:last-child,
.panel > .table-responsive:first-child > .table:first-child > thead:first-child > tr:first-child th:last-child,
.panel > .table:first-child > tbody:first-child > tr:first-child th:last-child,
.panel > .table-responsive:first-child > .table:first-child > tbody:first-child > tr:first-child th:last-child {
  border-top-right-radius: 1px;
}
.panel > .table:last-child,
.panel > .table-responsive:last-child > .table:last-child {
  border-bottom-right-radius: 1px;
  border-bottom-left-radius: 1px;
}
.panel > .table:last-child > tbody:last-child > tr:last-child,
.panel > .table-responsive:last-child > .table:last-child > tbody:last-child > tr:last-child,
.panel > .table:last-child > tfoot:last-child > tr:last-child,
.panel > .table-responsive:last-child > .table:last-child > tfoot:last-child > tr:last-child {
  border-bottom-left-radius: 1px;
  border-bottom-right-radius: 1px;
}
.panel > .table:last-child > tbody:last-child > tr:last-child td:first-child,
.panel > .table-responsive:last-child > .table:last-child > tbody:last-child > tr:last-child td:first-child,
.panel > .table:last-child > tfoot:last-child > tr:last-child td:first-child,
.panel > .table-responsive:last-child > .table:last-child > tfoot:last-child > tr:last-child td:first-child,
.panel > .table:last-child > tbody:last-child > tr:last-child th:first-child,
.panel > .table-responsive:last-child > .table:last-child > tbody:last-child > tr:last-child th:first-child,
.panel > .table:last-child > tfoot:last-child > tr:last-child th:first-child,
.panel > .table-responsive:last-child > .table:last-child > tfoot:last-child > tr:last-child th:first-child {
  border-bottom-left-radius: 1px;
}
.panel > .table:last-child > tbody:last-child > tr:last-child td:last-child,
.panel > .table-responsive:last-child > .table:last-child > tbody:last-child > tr:last-child td:last-child,
.panel > .table:last-child > tfoot:last-child > tr:last-child td:last-child,
.panel > .table-responsive:last-child > .table:last-child > tfoot:last-child > tr:last-child td:last-child,
.panel > .table:last-child > tbody:last-child > tr:last-child th:last-child,
.panel > .table-responsive:last-child > .table:last-child > tbody:last-child > tr:last-child th:last-child,
.panel > .table:last-child > tfoot:last-child > tr:last-child th:last-child,
.panel > .table-responsive:last-child > .table:last-child > tfoot:last-child > tr:last-child th:last-child {
  border-bottom-right-radius: 1px;
}
.panel > .panel-body + .table,
.panel > .panel-body + .table-responsive,
.panel > .table + .panel-body,
.panel > .table-responsive + .panel-body {
  border-top: 1px solid #ddd;
}
.panel > .table > tbody:first-child > tr:first-child th,
.panel > .table > tbody:first-child > tr:first-child td {
  border-top: 0;
}
.panel > .table-bordered,
.panel > .table-responsive > .table-bordered {
  border: 0;
}
.panel > .table-bordered > thead > tr > th:first-child,
.panel > .table-responsive > .table-bordered > thead > tr > th:first-child,
.panel > .table-bordered > tbody > tr > th:first-child,
.panel > .table-responsive > .table-bordered > tbody > tr > th:first-child,
.panel > .table-bordered > tfoot > tr > th:first-child,
.panel > .table-responsive > .table-bordered > tfoot > tr > th:first-child,
.panel > .table-bordered > thead > tr > td:first-child,
.panel > .table-responsive > .table-bordered > thead > tr > td:first-child,
.panel > .table-bordered > tbody > tr > td:first-child,
.panel > .table-responsive > .table-bordered > tbody > tr > td:first-child,
.panel > .table-bordered > tfoot > tr > td:first-child,
.panel > .table-responsive > .table-bordered > tfoot > tr > td:first-child {
  border-left: 0;
}
.panel > .table-bordered > thead > tr > th:last-child,
.panel > .table-responsive > .table-bordered > thead > tr > th:last-child,
.panel > .table-bordered > tbody > tr > th:last-child,
.panel > .table-responsive > .table-bordered > tbody > tr > th:last-child,
.panel > .table-bordered > tfoot > tr > th:last-child,
.panel > .table-responsive > .table-bordered > tfoot > tr > th:last-child,
.panel > .table-bordered > thead > tr > td:last-child,
.panel > .table-responsive > .table-bordered > thead > tr > td:last-child,
.panel > .table-bordered > tbody > tr > td:last-child,
.panel > .table-responsive > .table-bordered > tbody > tr > td:last-child,
.panel > .table-bordered > tfoot > tr > td:last-child,
.panel > .table-responsive > .table-bordered > tfoot > tr > td:last-child {
  border-right: 0;
}
.panel > .table-bordered > thead > tr:first-child > td,
.panel > .table-responsive > .table-bordered > thead > tr:first-child > td,
.panel > .table-bordered > tbody > tr:first-child > td,
.panel > .table-responsive > .table-bordered > tbody > tr:first-child > td,
.panel > .table-bordered > thead > tr:first-child > th,
.panel > .table-responsive > .table-bordered > thead > tr:first-child > th,
.panel > .table-bordered > tbody > tr:first-child > th,
.panel > .table-responsive > .table-bordered > tbody > tr:first-child > th {
  border-bottom: 0;
}
.panel > .table-bordered > tbody > tr:last-child > td,
.panel > .table-responsive > .table-bordered > tbody > tr:last-child > td,
.panel > .table-bordered > tfoot > tr:last-child > td,
.panel > .table-responsive > .table-bordered > tfoot > tr:last-child > td,
.panel > .table-bordered > tbody > tr:last-child > th,
.panel > .table-responsive > .table-bordered > tbody > tr:last-child > th,
.panel > .table-bordered > tfoot > tr:last-child > th,
.panel > .table-responsive > .table-bordered > tfoot > tr:last-child > th {
  border-bottom: 0;
}
.panel > .table-responsive {
  border: 0;
  margin-bottom: 0;
}
.panel-group {
  margin-bottom: 18px;
}
.panel-group .panel {
  margin-bottom: 0;
  border-radius: 2px;
}
.panel-group .panel + .panel {
  margin-top: 5px;
}
.panel-group .panel-heading {
  border-bottom: 0;
}
.panel-group .panel-heading + .panel-collapse > .panel-body,
.panel-group .panel-heading + .panel-collapse > .list-group {
  border-top: 1px solid #ddd;
}
.panel-group .panel-footer {
  border-top: 0;
}
.panel-group .panel-footer + .panel-collapse .panel-body {
  border-bottom: 1px solid #ddd;
}
.panel-default {
  border-color: #ddd;
}
.panel-default > .panel-heading {
  color: #333333;
  background-color: #f5f5f5;
  border-color: #ddd;
}
.panel-default > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #ddd;
}
.panel-default > .panel-heading .badge {
  color: #f5f5f5;
  background-color: #333333;
}
.panel-default > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #ddd;
}
.panel-primary {
  border-color: #337ab7;
}
.panel-primary > .panel-heading {
  color: #fff;
  background-color: #337ab7;
  border-color: #337ab7;
}
.panel-primary > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #337ab7;
}
.panel-primary > .panel-heading .badge {
  color: #337ab7;
  background-color: #fff;
}
.panel-primary > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #337ab7;
}
.panel-success {
  border-color: #d6e9c6;
}
.panel-success > .panel-heading {
  color: #3c763d;
  background-color: #dff0d8;
  border-color: #d6e9c6;
}
.panel-success > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #d6e9c6;
}
.panel-success > .panel-heading .badge {
  color: #dff0d8;
  background-color: #3c763d;
}
.panel-success > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #d6e9c6;
}
.panel-info {
  border-color: #bce8f1;
}
.panel-info > .panel-heading {
  color: #31708f;
  background-color: #d9edf7;
  border-color: #bce8f1;
}
.panel-info > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #bce8f1;
}
.panel-info > .panel-heading .badge {
  color: #d9edf7;
  background-color: #31708f;
}
.panel-info > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #bce8f1;
}
.panel-warning {
  border-color: #faebcc;
}
.panel-warning > .panel-heading {
  color: #8a6d3b;
  background-color: #fcf8e3;
  border-color: #faebcc;
}
.panel-warning > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #faebcc;
}
.panel-warning > .panel-heading .badge {
  color: #fcf8e3;
  background-color: #8a6d3b;
}
.panel-warning > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #faebcc;
}
.panel-danger {
  border-color: #ebccd1;
}
.panel-danger > .panel-heading {
  color: #a94442;
  background-color: #f2dede;
  border-color: #ebccd1;
}
.panel-danger > .panel-heading + .panel-collapse > .panel-body {
  border-top-color: #ebccd1;
}
.panel-danger > .panel-heading .badge {
  color: #f2dede;
  background-color: #a94442;
}
.panel-danger > .panel-footer + .panel-collapse > .panel-body {
  border-bottom-color: #ebccd1;
}
.embed-responsive {
  position: relative;
  display: block;
  height: 0;
  padding: 0;
  overflow: hidden;
}
.embed-responsive .embed-responsive-item,
.embed-responsive iframe,
.embed-responsive embed,
.embed-responsive object,
.embed-responsive video {
  position: absolute;
  top: 0;
  left: 0;
  bottom: 0;
  height: 100%;
  width: 100%;
  border: 0;
}
.embed-responsive-16by9 {
  padding-bottom: 56.25%;
}
.embed-responsive-4by3 {
  padding-bottom: 75%;
}
.well {
  min-height: 20px;
  padding: 19px;
  margin-bottom: 20px;
  background-color: #f5f5f5;
  border: 1px solid #e3e3e3;
  border-radius: 2px;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.05);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.05);
}
.well blockquote {
  border-color: #ddd;
  border-color: rgba(0, 0, 0, 0.15);
}
.well-lg {
  padding: 24px;
  border-radius: 3px;
}
.well-sm {
  padding: 9px;
  border-radius: 1px;
}
.close {
  float: right;
  font-size: 19.5px;
  font-weight: bold;
  line-height: 1;
  color: #000;
  text-shadow: 0 1px 0 #fff;
  opacity: 0.2;
  filter: alpha(opacity=20);
}
.close:hover,
.close:focus {
  color: #000;
  text-decoration: none;
  cursor: pointer;
  opacity: 0.5;
  filter: alpha(opacity=50);
}
button.close {
  padding: 0;
  cursor: pointer;
  background: transparent;
  border: 0;
  -webkit-appearance: none;
}
.modal-open {
  overflow: hidden;
}
.modal {
  display: none;
  overflow: hidden;
  position: fixed;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  z-index: 1050;
  -webkit-overflow-scrolling: touch;
  outline: 0;
}
.modal.fade .modal-dialog {
  -webkit-transform: translate(0, -25%);
  -ms-transform: translate(0, -25%);
  -o-transform: translate(0, -25%);
  transform: translate(0, -25%);
  -webkit-transition: -webkit-transform 0.3s ease-out;
  -moz-transition: -moz-transform 0.3s ease-out;
  -o-transition: -o-transform 0.3s ease-out;
  transition: transform 0.3s ease-out;
}
.modal.in .modal-dialog {
  -webkit-transform: translate(0, 0);
  -ms-transform: translate(0, 0);
  -o-transform: translate(0, 0);
  transform: translate(0, 0);
}
.modal-open .modal {
  overflow-x: hidden;
  overflow-y: auto;
}
.modal-dialog {
  position: relative;
  width: auto;
  margin: 10px;
}
.modal-content {
  position: relative;
  background-color: #fff;
  border: 1px solid #999;
  border: 1px solid rgba(0, 0, 0, 0.2);
  border-radius: 3px;
  -webkit-box-shadow: 0 3px 9px rgba(0, 0, 0, 0.5);
  box-shadow: 0 3px 9px rgba(0, 0, 0, 0.5);
  background-clip: padding-box;
  outline: 0;
}
.modal-backdrop {
  position: fixed;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  z-index: 1040;
  background-color: #000;
}
.modal-backdrop.fade {
  opacity: 0;
  filter: alpha(opacity=0);
}
.modal-backdrop.in {
  opacity: 0.5;
  filter: alpha(opacity=50);
}
.modal-header {
  padding: 15px;
  border-bottom: 1px solid #e5e5e5;
}
.modal-header .close {
  margin-top: -2px;
}
.modal-title {
  margin: 0;
  line-height: 1.42857143;
}
.modal-body {
  position: relative;
  padding: 15px;
}
.modal-footer {
  padding: 15px;
  text-align: right;
  border-top: 1px solid #e5e5e5;
}
.modal-footer .btn + .btn {
  margin-left: 5px;
  margin-bottom: 0;
}
.modal-footer .btn-group .btn + .btn {
  margin-left: -1px;
}
.modal-footer .btn-block + .btn-block {
  margin-left: 0;
}
.modal-scrollbar-measure {
  position: absolute;
  top: -9999px;
  width: 50px;
  height: 50px;
  overflow: scroll;
}
@media (min-width: 768px) {
  .modal-dialog {
    width: 600px;
    margin: 30px auto;
  }
  .modal-content {
    -webkit-box-shadow: 0 5px 15px rgba(0, 0, 0, 0.5);
    box-shadow: 0 5px 15px rgba(0, 0, 0, 0.5);
  }
  .modal-sm {
    width: 300px;
  }
}
@media (min-width: 992px) {
  .modal-lg {
    width: 900px;
  }
}
.tooltip {
  position: absolute;
  z-index: 1070;
  display: block;
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  font-style: normal;
  font-weight: normal;
  letter-spacing: normal;
  line-break: auto;
  line-height: 1.42857143;
  text-align: left;
  text-align: start;
  text-decoration: none;
  text-shadow: none;
  text-transform: none;
  white-space: normal;
  word-break: normal;
  word-spacing: normal;
  word-wrap: normal;
  font-size: 12px;
  opacity: 0;
  filter: alpha(opacity=0);
}
.tooltip.in {
  opacity: 0.9;
  filter: alpha(opacity=90);
}
.tooltip.top {
  margin-top: -3px;
  padding: 5px 0;
}
.tooltip.right {
  margin-left: 3px;
  padding: 0 5px;
}
.tooltip.bottom {
  margin-top: 3px;
  padding: 5px 0;
}
.tooltip.left {
  margin-left: -3px;
  padding: 0 5px;
}
.tooltip-inner {
  max-width: 200px;
  padding: 3px 8px;
  color: #fff;
  text-align: center;
  background-color: #000;
  border-radius: 2px;
}
.tooltip-arrow {
  position: absolute;
  width: 0;
  height: 0;
  border-color: transparent;
  border-style: solid;
}
.tooltip.top .tooltip-arrow {
  bottom: 0;
  left: 50%;
  margin-left: -5px;
  border-width: 5px 5px 0;
  border-top-color: #000;
}
.tooltip.top-left .tooltip-arrow {
  bottom: 0;
  right: 5px;
  margin-bottom: -5px;
  border-width: 5px 5px 0;
  border-top-color: #000;
}
.tooltip.top-right .tooltip-arrow {
  bottom: 0;
  left: 5px;
  margin-bottom: -5px;
  border-width: 5px 5px 0;
  border-top-color: #000;
}
.tooltip.right .tooltip-arrow {
  top: 50%;
  left: 0;
  margin-top: -5px;
  border-width: 5px 5px 5px 0;
  border-right-color: #000;
}
.tooltip.left .tooltip-arrow {
  top: 50%;
  right: 0;
  margin-top: -5px;
  border-width: 5px 0 5px 5px;
  border-left-color: #000;
}
.tooltip.bottom .tooltip-arrow {
  top: 0;
  left: 50%;
  margin-left: -5px;
  border-width: 0 5px 5px;
  border-bottom-color: #000;
}
.tooltip.bottom-left .tooltip-arrow {
  top: 0;
  right: 5px;
  margin-top: -5px;
  border-width: 0 5px 5px;
  border-bottom-color: #000;
}
.tooltip.bottom-right .tooltip-arrow {
  top: 0;
  left: 5px;
  margin-top: -5px;
  border-width: 0 5px 5px;
  border-bottom-color: #000;
}
.popover {
  position: absolute;
  top: 0;
  left: 0;
  z-index: 1060;
  display: none;
  max-width: 276px;
  padding: 1px;
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  font-style: normal;
  font-weight: normal;
  letter-spacing: normal;
  line-break: auto;
  line-height: 1.42857143;
  text-align: left;
  text-align: start;
  text-decoration: none;
  text-shadow: none;
  text-transform: none;
  white-space: normal;
  word-break: normal;
  word-spacing: normal;
  word-wrap: normal;
  font-size: 13px;
  background-color: #fff;
  background-clip: padding-box;
  border: 1px solid #ccc;
  border: 1px solid rgba(0, 0, 0, 0.2);
  border-radius: 3px;
  -webkit-box-shadow: 0 5px 10px rgba(0, 0, 0, 0.2);
  box-shadow: 0 5px 10px rgba(0, 0, 0, 0.2);
}
.popover.top {
  margin-top: -10px;
}
.popover.right {
  margin-left: 10px;
}
.popover.bottom {
  margin-top: 10px;
}
.popover.left {
  margin-left: -10px;
}
.popover-title {
  margin: 0;
  padding: 8px 14px;
  font-size: 13px;
  background-color: #f7f7f7;
  border-bottom: 1px solid #ebebeb;
  border-radius: 2px 2px 0 0;
}
.popover-content {
  padding: 9px 14px;
}
.popover > .arrow,
.popover > .arrow:after {
  position: absolute;
  display: block;
  width: 0;
  height: 0;
  border-color: transparent;
  border-style: solid;
}
.popover > .arrow {
  border-width: 11px;
}
.popover > .arrow:after {
  border-width: 10px;
  content: "";
}
.popover.top > .arrow {
  left: 50%;
  margin-left: -11px;
  border-bottom-width: 0;
  border-top-color: #999999;
  border-top-color: rgba(0, 0, 0, 0.25);
  bottom: -11px;
}
.popover.top > .arrow:after {
  content: " ";
  bottom: 1px;
  margin-left: -10px;
  border-bottom-width: 0;
  border-top-color: #fff;
}
.popover.right > .arrow {
  top: 50%;
  left: -11px;
  margin-top: -11px;
  border-left-width: 0;
  border-right-color: #999999;
  border-right-color: rgba(0, 0, 0, 0.25);
}
.popover.right > .arrow:after {
  content: " ";
  left: 1px;
  bottom: -10px;
  border-left-width: 0;
  border-right-color: #fff;
}
.popover.bottom > .arrow {
  left: 50%;
  margin-left: -11px;
  border-top-width: 0;
  border-bottom-color: #999999;
  border-bottom-color: rgba(0, 0, 0, 0.25);
  top: -11px;
}
.popover.bottom > .arrow:after {
  content: " ";
  top: 1px;
  margin-left: -10px;
  border-top-width: 0;
  border-bottom-color: #fff;
}
.popover.left > .arrow {
  top: 50%;
  right: -11px;
  margin-top: -11px;
  border-right-width: 0;
  border-left-color: #999999;
  border-left-color: rgba(0, 0, 0, 0.25);
}
.popover.left > .arrow:after {
  content: " ";
  right: 1px;
  border-right-width: 0;
  border-left-color: #fff;
  bottom: -10px;
}
.carousel {
  position: relative;
}
.carousel-inner {
  position: relative;
  overflow: hidden;
  width: 100%;
}
.carousel-inner > .item {
  display: none;
  position: relative;
  -webkit-transition: 0.6s ease-in-out left;
  -o-transition: 0.6s ease-in-out left;
  transition: 0.6s ease-in-out left;
}
.carousel-inner > .item > img,
.carousel-inner > .item > a > img {
  line-height: 1;
}
@media all and (transform-3d), (-webkit-transform-3d) {
  .carousel-inner > .item {
    -webkit-transition: -webkit-transform 0.6s ease-in-out;
    -moz-transition: -moz-transform 0.6s ease-in-out;
    -o-transition: -o-transform 0.6s ease-in-out;
    transition: transform 0.6s ease-in-out;
    -webkit-backface-visibility: hidden;
    -moz-backface-visibility: hidden;
    backface-visibility: hidden;
    -webkit-perspective: 1000px;
    -moz-perspective: 1000px;
    perspective: 1000px;
  }
  .carousel-inner > .item.next,
  .carousel-inner > .item.active.right {
    -webkit-transform: translate3d(100%, 0, 0);
    transform: translate3d(100%, 0, 0);
    left: 0;
  }
  .carousel-inner > .item.prev,
  .carousel-inner > .item.active.left {
    -webkit-transform: translate3d(-100%, 0, 0);
    transform: translate3d(-100%, 0, 0);
    left: 0;
  }
  .carousel-inner > .item.next.left,
  .carousel-inner > .item.prev.right,
  .carousel-inner > .item.active {
    -webkit-transform: translate3d(0, 0, 0);
    transform: translate3d(0, 0, 0);
    left: 0;
  }
}
.carousel-inner > .active,
.carousel-inner > .next,
.carousel-inner > .prev {
  display: block;
}
.carousel-inner > .active {
  left: 0;
}
.carousel-inner > .next,
.carousel-inner > .prev {
  position: absolute;
  top: 0;
  width: 100%;
}
.carousel-inner > .next {
  left: 100%;
}
.carousel-inner > .prev {
  left: -100%;
}
.carousel-inner > .next.left,
.carousel-inner > .prev.right {
  left: 0;
}
.carousel-inner > .active.left {
  left: -100%;
}
.carousel-inner > .active.right {
  left: 100%;
}
.carousel-control {
  position: absolute;
  top: 0;
  left: 0;
  bottom: 0;
  width: 15%;
  opacity: 0.5;
  filter: alpha(opacity=50);
  font-size: 20px;
  color: #fff;
  text-align: center;
  text-shadow: 0 1px 2px rgba(0, 0, 0, 0.6);
  background-color: rgba(0, 0, 0, 0);
}
.carousel-control.left {
  background-image: -webkit-linear-gradient(left, rgba(0, 0, 0, 0.5) 0%, rgba(0, 0, 0, 0.0001) 100%);
  background-image: -o-linear-gradient(left, rgba(0, 0, 0, 0.5) 0%, rgba(0, 0, 0, 0.0001) 100%);
  background-image: linear-gradient(to right, rgba(0, 0, 0, 0.5) 0%, rgba(0, 0, 0, 0.0001) 100%);
  background-repeat: repeat-x;
  filter: progid:DXImageTransform.Microsoft.gradient(startColorstr='#80000000', endColorstr='#00000000', GradientType=1);
}
.carousel-control.right {
  left: auto;
  right: 0;
  background-image: -webkit-linear-gradient(left, rgba(0, 0, 0, 0.0001) 0%, rgba(0, 0, 0, 0.5) 100%);
  background-image: -o-linear-gradient(left, rgba(0, 0, 0, 0.0001) 0%, rgba(0, 0, 0, 0.5) 100%);
  background-image: linear-gradient(to right, rgba(0, 0, 0, 0.0001) 0%, rgba(0, 0, 0, 0.5) 100%);
  background-repeat: repeat-x;
  filter: progid:DXImageTransform.Microsoft.gradient(startColorstr='#00000000', endColorstr='#80000000', GradientType=1);
}
.carousel-control:hover,
.carousel-control:focus {
  outline: 0;
  color: #fff;
  text-decoration: none;
  opacity: 0.9;
  filter: alpha(opacity=90);
}
.carousel-control .icon-prev,
.carousel-control .icon-next,
.carousel-control .glyphicon-chevron-left,
.carousel-control .glyphicon-chevron-right {
  position: absolute;
  top: 50%;
  margin-top: -10px;
  z-index: 5;
  display: inline-block;
}
.carousel-control .icon-prev,
.carousel-control .glyphicon-chevron-left {
  left: 50%;
  margin-left: -10px;
}
.carousel-control .icon-next,
.carousel-control .glyphicon-chevron-right {
  right: 50%;
  margin-right: -10px;
}
.carousel-control .icon-prev,
.carousel-control .icon-next {
  width: 20px;
  height: 20px;
  line-height: 1;
  font-family: serif;
}
.carousel-control .icon-prev:before {
  content: '\2039';
}
.carousel-control .icon-next:before {
  content: '\203a';
}
.carousel-indicators {
  position: absolute;
  bottom: 10px;
  left: 50%;
  z-index: 15;
  width: 60%;
  margin-left: -30%;
  padding-left: 0;
  list-style: none;
  text-align: center;
}
.carousel-indicators li {
  display: inline-block;
  width: 10px;
  height: 10px;
  margin: 1px;
  text-indent: -999px;
  border: 1px solid #fff;
  border-radius: 10px;
  cursor: pointer;
  background-color: #000 \9;
  background-color: rgba(0, 0, 0, 0);
}
.carousel-indicators .active {
  margin: 0;
  width: 12px;
  height: 12px;
  background-color: #fff;
}
.carousel-caption {
  position: absolute;
  left: 15%;
  right: 15%;
  bottom: 20px;
  z-index: 10;
  padding-top: 20px;
  padding-bottom: 20px;
  color: #fff;
  text-align: center;
  text-shadow: 0 1px 2px rgba(0, 0, 0, 0.6);
}
.carousel-caption .btn {
  text-shadow: none;
}
@media screen and (min-width: 768px) {
  .carousel-control .glyphicon-chevron-left,
  .carousel-control .glyphicon-chevron-right,
  .carousel-control .icon-prev,
  .carousel-control .icon-next {
    width: 30px;
    height: 30px;
    margin-top: -10px;
    font-size: 30px;
  }
  .carousel-control .glyphicon-chevron-left,
  .carousel-control .icon-prev {
    margin-left: -10px;
  }
  .carousel-control .glyphicon-chevron-right,
  .carousel-control .icon-next {
    margin-right: -10px;
  }
  .carousel-caption {
    left: 20%;
    right: 20%;
    padding-bottom: 30px;
  }
  .carousel-indicators {
    bottom: 20px;
  }
}
.clearfix:before,
.clearfix:after,
.dl-horizontal dd:before,
.dl-horizontal dd:after,
.container:before,
.container:after,
.container-fluid:before,
.container-fluid:after,
.row:before,
.row:after,
.form-horizontal .form-group:before,
.form-horizontal .form-group:after,
.btn-toolbar:before,
.btn-toolbar:after,
.btn-group-vertical > .btn-group:before,
.btn-group-vertical > .btn-group:after,
.nav:before,
.nav:after,
.navbar:before,
.navbar:after,
.navbar-header:before,
.navbar-header:after,
.navbar-collapse:before,
.navbar-collapse:after,
.pager:before,
.pager:after,
.panel-body:before,
.panel-body:after,
.modal-header:before,
.modal-header:after,
.modal-footer:before,
.modal-footer:after,
.item_buttons:before,
.item_buttons:after {
  content: " ";
  display: table;
}
.clearfix:after,
.dl-horizontal dd:after,
.container:after,
.container-fluid:after,
.row:after,
.form-horizontal .form-group:after,
.btn-toolbar:after,
.btn-group-vertical > .btn-group:after,
.nav:after,
.navbar:after,
.navbar-header:after,
.navbar-collapse:after,
.pager:after,
.panel-body:after,
.modal-header:after,
.modal-footer:after,
.item_buttons:after {
  clear: both;
}
.center-block {
  display: block;
  margin-left: auto;
  margin-right: auto;
}
.pull-right {
  float: right !important;
}
.pull-left {
  float: left !important;
}
.hide {
  display: none !important;
}
.show {
  display: block !important;
}
.invisible {
  visibility: hidden;
}
.text-hide {
  font: 0/0 a;
  color: transparent;
  text-shadow: none;
  background-color: transparent;
  border: 0;
}
.hidden {
  display: none !important;
}
.affix {
  position: fixed;
}
@-ms-viewport {
  width: device-width;
}
.visible-xs,
.visible-sm,
.visible-md,
.visible-lg {
  display: none !important;
}
.visible-xs-block,
.visible-xs-inline,
.visible-xs-inline-block,
.visible-sm-block,
.visible-sm-inline,
.visible-sm-inline-block,
.visible-md-block,
.visible-md-inline,
.visible-md-inline-block,
.visible-lg-block,
.visible-lg-inline,
.visible-lg-inline-block {
  display: none !important;
}
@media (max-width: 767px) {
  .visible-xs {
    display: block !important;
  }
  table.visible-xs {
    display: table !important;
  }
  tr.visible-xs {
    display: table-row !important;
  }
  th.visible-xs,
  td.visible-xs {
    display: table-cell !important;
  }
}
@media (max-width: 767px) {
  .visible-xs-block {
    display: block !important;
  }
}
@media (max-width: 767px) {
  .visible-xs-inline {
    display: inline !important;
  }
}
@media (max-width: 767px) {
  .visible-xs-inline-block {
    display: inline-block !important;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  .visible-sm {
    display: block !important;
  }
  table.visible-sm {
    display: table !important;
  }
  tr.visible-sm {
    display: table-row !important;
  }
  th.visible-sm,
  td.visible-sm {
    display: table-cell !important;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  .visible-sm-block {
    display: block !important;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  .visible-sm-inline {
    display: inline !important;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  .visible-sm-inline-block {
    display: inline-block !important;
  }
}
@media (min-width: 992px) and (max-width: 1199px) {
  .visible-md {
    display: block !important;
  }
  table.visible-md {
    display: table !important;
  }
  tr.visible-md {
    display: table-row !important;
  }
  th.visible-md,
  td.visible-md {
    display: table-cell !important;
  }
}
@media (min-width: 992px) and (max-width: 1199px) {
  .visible-md-block {
    display: block !important;
  }
}
@media (min-width: 992px) and (max-width: 1199px) {
  .visible-md-inline {
    display: inline !important;
  }
}
@media (min-width: 992px) and (max-width: 1199px) {
  .visible-md-inline-block {
    display: inline-block !important;
  }
}
@media (min-width: 1200px) {
  .visible-lg {
    display: block !important;
  }
  table.visible-lg {
    display: table !important;
  }
  tr.visible-lg {
    display: table-row !important;
  }
  th.visible-lg,
  td.visible-lg {
    display: table-cell !important;
  }
}
@media (min-width: 1200px) {
  .visible-lg-block {
    display: block !important;
  }
}
@media (min-width: 1200px) {
  .visible-lg-inline {
    display: inline !important;
  }
}
@media (min-width: 1200px) {
  .visible-lg-inline-block {
    display: inline-block !important;
  }
}
@media (max-width: 767px) {
  .hidden-xs {
    display: none !important;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  .hidden-sm {
    display: none !important;
  }
}
@media (min-width: 992px) and (max-width: 1199px) {
  .hidden-md {
    display: none !important;
  }
}
@media (min-width: 1200px) {
  .hidden-lg {
    display: none !important;
  }
}
.visible-print {
  display: none !important;
}
@media print {
  .visible-print {
    display: block !important;
  }
  table.visible-print {
    display: table !important;
  }
  tr.visible-print {
    display: table-row !important;
  }
  th.visible-print,
  td.visible-print {
    display: table-cell !important;
  }
}
.visible-print-block {
  display: none !important;
}
@media print {
  .visible-print-block {
    display: block !important;
  }
}
.visible-print-inline {
  display: none !important;
}
@media print {
  .visible-print-inline {
    display: inline !important;
  }
}
.visible-print-inline-block {
  display: none !important;
}
@media print {
  .visible-print-inline-block {
    display: inline-block !important;
  }
}
@media print {
  .hidden-print {
    display: none !important;
  }
}
/*!
*
* Font Awesome
*
*/
/*!
 *  Font Awesome 4.2.0 by @davegandy - http://fontawesome.io - @fontawesome
 *  License - http://fontawesome.io/license (Font: SIL OFL 1.1, CSS: MIT License)
 */
/* FONT PATH
 * -------------------------- */
@font-face {
  font-family: 'FontAwesome';
  src: url('../components/font-awesome/fonts/fontawesome-webfont.eot?v=4.2.0');
  src: url('../components/font-awesome/fonts/fontawesome-webfont.eot?#iefix&v=4.2.0') format('embedded-opentype'), url('../components/font-awesome/fonts/fontawesome-webfont.woff?v=4.2.0') format('woff'), url('../components/font-awesome/fonts/fontawesome-webfont.ttf?v=4.2.0') format('truetype'), url('../components/font-awesome/fonts/fontawesome-webfont.svg?v=4.2.0#fontawesomeregular') format('svg');
  font-weight: normal;
  font-style: normal;
}
.fa {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
/* makes the font 33% larger relative to the icon container */
.fa-lg {
  font-size: 1.33333333em;
  line-height: 0.75em;
  vertical-align: -15%;
}
.fa-2x {
  font-size: 2em;
}
.fa-3x {
  font-size: 3em;
}
.fa-4x {
  font-size: 4em;
}
.fa-5x {
  font-size: 5em;
}
.fa-fw {
  width: 1.28571429em;
  text-align: center;
}
.fa-ul {
  padding-left: 0;
  margin-left: 2.14285714em;
  list-style-type: none;
}
.fa-ul > li {
  position: relative;
}
.fa-li {
  position: absolute;
  left: -2.14285714em;
  width: 2.14285714em;
  top: 0.14285714em;
  text-align: center;
}
.fa-li.fa-lg {
  left: -1.85714286em;
}
.fa-border {
  padding: .2em .25em .15em;
  border: solid 0.08em #eee;
  border-radius: .1em;
}
.pull-right {
  float: right;
}
.pull-left {
  float: left;
}
.fa.pull-left {
  margin-right: .3em;
}
.fa.pull-right {
  margin-left: .3em;
}
.fa-spin {
  -webkit-animation: fa-spin 2s infinite linear;
  animation: fa-spin 2s infinite linear;
}
@-webkit-keyframes fa-spin {
  0% {
    -webkit-transform: rotate(0deg);
    transform: rotate(0deg);
  }
  100% {
    -webkit-transform: rotate(359deg);
    transform: rotate(359deg);
  }
}
@keyframes fa-spin {
  0% {
    -webkit-transform: rotate(0deg);
    transform: rotate(0deg);
  }
  100% {
    -webkit-transform: rotate(359deg);
    transform: rotate(359deg);
  }
}
.fa-rotate-90 {
  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=1);
  -webkit-transform: rotate(90deg);
  -ms-transform: rotate(90deg);
  transform: rotate(90deg);
}
.fa-rotate-180 {
  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=2);
  -webkit-transform: rotate(180deg);
  -ms-transform: rotate(180deg);
  transform: rotate(180deg);
}
.fa-rotate-270 {
  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=3);
  -webkit-transform: rotate(270deg);
  -ms-transform: rotate(270deg);
  transform: rotate(270deg);
}
.fa-flip-horizontal {
  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=0, mirror=1);
  -webkit-transform: scale(-1, 1);
  -ms-transform: scale(-1, 1);
  transform: scale(-1, 1);
}
.fa-flip-vertical {
  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=2, mirror=1);
  -webkit-transform: scale(1, -1);
  -ms-transform: scale(1, -1);
  transform: scale(1, -1);
}
:root .fa-rotate-90,
:root .fa-rotate-180,
:root .fa-rotate-270,
:root .fa-flip-horizontal,
:root .fa-flip-vertical {
  filter: none;
}
.fa-stack {
  position: relative;
  display: inline-block;
  width: 2em;
  height: 2em;
  line-height: 2em;
  vertical-align: middle;
}
.fa-stack-1x,
.fa-stack-2x {
  position: absolute;
  left: 0;
  width: 100%;
  text-align: center;
}
.fa-stack-1x {
  line-height: inherit;
}
.fa-stack-2x {
  font-size: 2em;
}
.fa-inverse {
  color: #fff;
}
/* Font Awesome uses the Unicode Private Use Area (PUA) to ensure screen
   readers do not read off random characters that represent icons */
.fa-glass:before {
  content: "\f000";
}
.fa-music:before {
  content: "\f001";
}
.fa-search:before {
  content: "\f002";
}
.fa-envelope-o:before {
  content: "\f003";
}
.fa-heart:before {
  content: "\f004";
}
.fa-star:before {
  content: "\f005";
}
.fa-star-o:before {
  content: "\f006";
}
.fa-user:before {
  content: "\f007";
}
.fa-film:before {
  content: "\f008";
}
.fa-th-large:before {
  content: "\f009";
}
.fa-th:before {
  content: "\f00a";
}
.fa-th-list:before {
  content: "\f00b";
}
.fa-check:before {
  content: "\f00c";
}
.fa-remove:before,
.fa-close:before,
.fa-times:before {
  content: "\f00d";
}
.fa-search-plus:before {
  content: "\f00e";
}
.fa-search-minus:before {
  content: "\f010";
}
.fa-power-off:before {
  content: "\f011";
}
.fa-signal:before {
  content: "\f012";
}
.fa-gear:before,
.fa-cog:before {
  content: "\f013";
}
.fa-trash-o:before {
  content: "\f014";
}
.fa-home:before {
  content: "\f015";
}
.fa-file-o:before {
  content: "\f016";
}
.fa-clock-o:before {
  content: "\f017";
}
.fa-road:before {
  content: "\f018";
}
.fa-download:before {
  content: "\f019";
}
.fa-arrow-circle-o-down:before {
  content: "\f01a";
}
.fa-arrow-circle-o-up:before {
  content: "\f01b";
}
.fa-inbox:before {
  content: "\f01c";
}
.fa-play-circle-o:before {
  content: "\f01d";
}
.fa-rotate-right:before,
.fa-repeat:before {
  content: "\f01e";
}
.fa-refresh:before {
  content: "\f021";
}
.fa-list-alt:before {
  content: "\f022";
}
.fa-lock:before {
  content: "\f023";
}
.fa-flag:before {
  content: "\f024";
}
.fa-headphones:before {
  content: "\f025";
}
.fa-volume-off:before {
  content: "\f026";
}
.fa-volume-down:before {
  content: "\f027";
}
.fa-volume-up:before {
  content: "\f028";
}
.fa-qrcode:before {
  content: "\f029";
}
.fa-barcode:before {
  content: "\f02a";
}
.fa-tag:before {
  content: "\f02b";
}
.fa-tags:before {
  content: "\f02c";
}
.fa-book:before {
  content: "\f02d";
}
.fa-bookmark:before {
  content: "\f02e";
}
.fa-print:before {
  content: "\f02f";
}
.fa-camera:before {
  content: "\f030";
}
.fa-font:before {
  content: "\f031";
}
.fa-bold:before {
  content: "\f032";
}
.fa-italic:before {
  content: "\f033";
}
.fa-text-height:before {
  content: "\f034";
}
.fa-text-width:before {
  content: "\f035";
}
.fa-align-left:before {
  content: "\f036";
}
.fa-align-center:before {
  content: "\f037";
}
.fa-align-right:before {
  content: "\f038";
}
.fa-align-justify:before {
  content: "\f039";
}
.fa-list:before {
  content: "\f03a";
}
.fa-dedent:before,
.fa-outdent:before {
  content: "\f03b";
}
.fa-indent:before {
  content: "\f03c";
}
.fa-video-camera:before {
  content: "\f03d";
}
.fa-photo:before,
.fa-image:before,
.fa-picture-o:before {
  content: "\f03e";
}
.fa-pencil:before {
  content: "\f040";
}
.fa-map-marker:before {
  content: "\f041";
}
.fa-adjust:before {
  content: "\f042";
}
.fa-tint:before {
  content: "\f043";
}
.fa-edit:before,
.fa-pencil-square-o:before {
  content: "\f044";
}
.fa-share-square-o:before {
  content: "\f045";
}
.fa-check-square-o:before {
  content: "\f046";
}
.fa-arrows:before {
  content: "\f047";
}
.fa-step-backward:before {
  content: "\f048";
}
.fa-fast-backward:before {
  content: "\f049";
}
.fa-backward:before {
  content: "\f04a";
}
.fa-play:before {
  content: "\f04b";
}
.fa-pause:before {
  content: "\f04c";
}
.fa-stop:before {
  content: "\f04d";
}
.fa-forward:before {
  content: "\f04e";
}
.fa-fast-forward:before {
  content: "\f050";
}
.fa-step-forward:before {
  content: "\f051";
}
.fa-eject:before {
  content: "\f052";
}
.fa-chevron-left:before {
  content: "\f053";
}
.fa-chevron-right:before {
  content: "\f054";
}
.fa-plus-circle:before {
  content: "\f055";
}
.fa-minus-circle:before {
  content: "\f056";
}
.fa-times-circle:before {
  content: "\f057";
}
.fa-check-circle:before {
  content: "\f058";
}
.fa-question-circle:before {
  content: "\f059";
}
.fa-info-circle:before {
  content: "\f05a";
}
.fa-crosshairs:before {
  content: "\f05b";
}
.fa-times-circle-o:before {
  content: "\f05c";
}
.fa-check-circle-o:before {
  content: "\f05d";
}
.fa-ban:before {
  content: "\f05e";
}
.fa-arrow-left:before {
  content: "\f060";
}
.fa-arrow-right:before {
  content: "\f061";
}
.fa-arrow-up:before {
  content: "\f062";
}
.fa-arrow-down:before {
  content: "\f063";
}
.fa-mail-forward:before,
.fa-share:before {
  content: "\f064";
}
.fa-expand:before {
  content: "\f065";
}
.fa-compress:before {
  content: "\f066";
}
.fa-plus:before {
  content: "\f067";
}
.fa-minus:before {
  content: "\f068";
}
.fa-asterisk:before {
  content: "\f069";
}
.fa-exclamation-circle:before {
  content: "\f06a";
}
.fa-gift:before {
  content: "\f06b";
}
.fa-leaf:before {
  content: "\f06c";
}
.fa-fire:before {
  content: "\f06d";
}
.fa-eye:before {
  content: "\f06e";
}
.fa-eye-slash:before {
  content: "\f070";
}
.fa-warning:before,
.fa-exclamation-triangle:before {
  content: "\f071";
}
.fa-plane:before {
  content: "\f072";
}
.fa-calendar:before {
  content: "\f073";
}
.fa-random:before {
  content: "\f074";
}
.fa-comment:before {
  content: "\f075";
}
.fa-magnet:before {
  content: "\f076";
}
.fa-chevron-up:before {
  content: "\f077";
}
.fa-chevron-down:before {
  content: "\f078";
}
.fa-retweet:before {
  content: "\f079";
}
.fa-shopping-cart:before {
  content: "\f07a";
}
.fa-folder:before {
  content: "\f07b";
}
.fa-folder-open:before {
  content: "\f07c";
}
.fa-arrows-v:before {
  content: "\f07d";
}
.fa-arrows-h:before {
  content: "\f07e";
}
.fa-bar-chart-o:before,
.fa-bar-chart:before {
  content: "\f080";
}
.fa-twitter-square:before {
  content: "\f081";
}
.fa-facebook-square:before {
  content: "\f082";
}
.fa-camera-retro:before {
  content: "\f083";
}
.fa-key:before {
  content: "\f084";
}
.fa-gears:before,
.fa-cogs:before {
  content: "\f085";
}
.fa-comments:before {
  content: "\f086";
}
.fa-thumbs-o-up:before {
  content: "\f087";
}
.fa-thumbs-o-down:before {
  content: "\f088";
}
.fa-star-half:before {
  content: "\f089";
}
.fa-heart-o:before {
  content: "\f08a";
}
.fa-sign-out:before {
  content: "\f08b";
}
.fa-linkedin-square:before {
  content: "\f08c";
}
.fa-thumb-tack:before {
  content: "\f08d";
}
.fa-external-link:before {
  content: "\f08e";
}
.fa-sign-in:before {
  content: "\f090";
}
.fa-trophy:before {
  content: "\f091";
}
.fa-github-square:before {
  content: "\f092";
}
.fa-upload:before {
  content: "\f093";
}
.fa-lemon-o:before {
  content: "\f094";
}
.fa-phone:before {
  content: "\f095";
}
.fa-square-o:before {
  content: "\f096";
}
.fa-bookmark-o:before {
  content: "\f097";
}
.fa-phone-square:before {
  content: "\f098";
}
.fa-twitter:before {
  content: "\f099";
}
.fa-facebook:before {
  content: "\f09a";
}
.fa-github:before {
  content: "\f09b";
}
.fa-unlock:before {
  content: "\f09c";
}
.fa-credit-card:before {
  content: "\f09d";
}
.fa-rss:before {
  content: "\f09e";
}
.fa-hdd-o:before {
  content: "\f0a0";
}
.fa-bullhorn:before {
  content: "\f0a1";
}
.fa-bell:before {
  content: "\f0f3";
}
.fa-certificate:before {
  content: "\f0a3";
}
.fa-hand-o-right:before {
  content: "\f0a4";
}
.fa-hand-o-left:before {
  content: "\f0a5";
}
.fa-hand-o-up:before {
  content: "\f0a6";
}
.fa-hand-o-down:before {
  content: "\f0a7";
}
.fa-arrow-circle-left:before {
  content: "\f0a8";
}
.fa-arrow-circle-right:before {
  content: "\f0a9";
}
.fa-arrow-circle-up:before {
  content: "\f0aa";
}
.fa-arrow-circle-down:before {
  content: "\f0ab";
}
.fa-globe:before {
  content: "\f0ac";
}
.fa-wrench:before {
  content: "\f0ad";
}
.fa-tasks:before {
  content: "\f0ae";
}
.fa-filter:before {
  content: "\f0b0";
}
.fa-briefcase:before {
  content: "\f0b1";
}
.fa-arrows-alt:before {
  content: "\f0b2";
}
.fa-group:before,
.fa-users:before {
  content: "\f0c0";
}
.fa-chain:before,
.fa-link:before {
  content: "\f0c1";
}
.fa-cloud:before {
  content: "\f0c2";
}
.fa-flask:before {
  content: "\f0c3";
}
.fa-cut:before,
.fa-scissors:before {
  content: "\f0c4";
}
.fa-copy:before,
.fa-files-o:before {
  content: "\f0c5";
}
.fa-paperclip:before {
  content: "\f0c6";
}
.fa-save:before,
.fa-floppy-o:before {
  content: "\f0c7";
}
.fa-square:before {
  content: "\f0c8";
}
.fa-navicon:before,
.fa-reorder:before,
.fa-bars:before {
  content: "\f0c9";
}
.fa-list-ul:before {
  content: "\f0ca";
}
.fa-list-ol:before {
  content: "\f0cb";
}
.fa-strikethrough:before {
  content: "\f0cc";
}
.fa-underline:before {
  content: "\f0cd";
}
.fa-table:before {
  content: "\f0ce";
}
.fa-magic:before {
  content: "\f0d0";
}
.fa-truck:before {
  content: "\f0d1";
}
.fa-pinterest:before {
  content: "\f0d2";
}
.fa-pinterest-square:before {
  content: "\f0d3";
}
.fa-google-plus-square:before {
  content: "\f0d4";
}
.fa-google-plus:before {
  content: "\f0d5";
}
.fa-money:before {
  content: "\f0d6";
}
.fa-caret-down:before {
  content: "\f0d7";
}
.fa-caret-up:before {
  content: "\f0d8";
}
.fa-caret-left:before {
  content: "\f0d9";
}
.fa-caret-right:before {
  content: "\f0da";
}
.fa-columns:before {
  content: "\f0db";
}
.fa-unsorted:before,
.fa-sort:before {
  content: "\f0dc";
}
.fa-sort-down:before,
.fa-sort-desc:before {
  content: "\f0dd";
}
.fa-sort-up:before,
.fa-sort-asc:before {
  content: "\f0de";
}
.fa-envelope:before {
  content: "\f0e0";
}
.fa-linkedin:before {
  content: "\f0e1";
}
.fa-rotate-left:before,
.fa-undo:before {
  content: "\f0e2";
}
.fa-legal:before,
.fa-gavel:before {
  content: "\f0e3";
}
.fa-dashboard:before,
.fa-tachometer:before {
  content: "\f0e4";
}
.fa-comment-o:before {
  content: "\f0e5";
}
.fa-comments-o:before {
  content: "\f0e6";
}
.fa-flash:before,
.fa-bolt:before {
  content: "\f0e7";
}
.fa-sitemap:before {
  content: "\f0e8";
}
.fa-umbrella:before {
  content: "\f0e9";
}
.fa-paste:before,
.fa-clipboard:before {
  content: "\f0ea";
}
.fa-lightbulb-o:before {
  content: "\f0eb";
}
.fa-exchange:before {
  content: "\f0ec";
}
.fa-cloud-download:before {
  content: "\f0ed";
}
.fa-cloud-upload:before {
  content: "\f0ee";
}
.fa-user-md:before {
  content: "\f0f0";
}
.fa-stethoscope:before {
  content: "\f0f1";
}
.fa-suitcase:before {
  content: "\f0f2";
}
.fa-bell-o:before {
  content: "\f0a2";
}
.fa-coffee:before {
  content: "\f0f4";
}
.fa-cutlery:before {
  content: "\f0f5";
}
.fa-file-text-o:before {
  content: "\f0f6";
}
.fa-building-o:before {
  content: "\f0f7";
}
.fa-hospital-o:before {
  content: "\f0f8";
}
.fa-ambulance:before {
  content: "\f0f9";
}
.fa-medkit:before {
  content: "\f0fa";
}
.fa-fighter-jet:before {
  content: "\f0fb";
}
.fa-beer:before {
  content: "\f0fc";
}
.fa-h-square:before {
  content: "\f0fd";
}
.fa-plus-square:before {
  content: "\f0fe";
}
.fa-angle-double-left:before {
  content: "\f100";
}
.fa-angle-double-right:before {
  content: "\f101";
}
.fa-angle-double-up:before {
  content: "\f102";
}
.fa-angle-double-down:before {
  content: "\f103";
}
.fa-angle-left:before {
  content: "\f104";
}
.fa-angle-right:before {
  content: "\f105";
}
.fa-angle-up:before {
  content: "\f106";
}
.fa-angle-down:before {
  content: "\f107";
}
.fa-desktop:before {
  content: "\f108";
}
.fa-laptop:before {
  content: "\f109";
}
.fa-tablet:before {
  content: "\f10a";
}
.fa-mobile-phone:before,
.fa-mobile:before {
  content: "\f10b";
}
.fa-circle-o:before {
  content: "\f10c";
}
.fa-quote-left:before {
  content: "\f10d";
}
.fa-quote-right:before {
  content: "\f10e";
}
.fa-spinner:before {
  content: "\f110";
}
.fa-circle:before {
  content: "\f111";
}
.fa-mail-reply:before,
.fa-reply:before {
  content: "\f112";
}
.fa-github-alt:before {
  content: "\f113";
}
.fa-folder-o:before {
  content: "\f114";
}
.fa-folder-open-o:before {
  content: "\f115";
}
.fa-smile-o:before {
  content: "\f118";
}
.fa-frown-o:before {
  content: "\f119";
}
.fa-meh-o:before {
  content: "\f11a";
}
.fa-gamepad:before {
  content: "\f11b";
}
.fa-keyboard-o:before {
  content: "\f11c";
}
.fa-flag-o:before {
  content: "\f11d";
}
.fa-flag-checkered:before {
  content: "\f11e";
}
.fa-terminal:before {
  content: "\f120";
}
.fa-code:before {
  content: "\f121";
}
.fa-mail-reply-all:before,
.fa-reply-all:before {
  content: "\f122";
}
.fa-star-half-empty:before,
.fa-star-half-full:before,
.fa-star-half-o:before {
  content: "\f123";
}
.fa-location-arrow:before {
  content: "\f124";
}
.fa-crop:before {
  content: "\f125";
}
.fa-code-fork:before {
  content: "\f126";
}
.fa-unlink:before,
.fa-chain-broken:before {
  content: "\f127";
}
.fa-question:before {
  content: "\f128";
}
.fa-info:before {
  content: "\f129";
}
.fa-exclamation:before {
  content: "\f12a";
}
.fa-superscript:before {
  content: "\f12b";
}
.fa-subscript:before {
  content: "\f12c";
}
.fa-eraser:before {
  content: "\f12d";
}
.fa-puzzle-piece:before {
  content: "\f12e";
}
.fa-microphone:before {
  content: "\f130";
}
.fa-microphone-slash:before {
  content: "\f131";
}
.fa-shield:before {
  content: "\f132";
}
.fa-calendar-o:before {
  content: "\f133";
}
.fa-fire-extinguisher:before {
  content: "\f134";
}
.fa-rocket:before {
  content: "\f135";
}
.fa-maxcdn:before {
  content: "\f136";
}
.fa-chevron-circle-left:before {
  content: "\f137";
}
.fa-chevron-circle-right:before {
  content: "\f138";
}
.fa-chevron-circle-up:before {
  content: "\f139";
}
.fa-chevron-circle-down:before {
  content: "\f13a";
}
.fa-html5:before {
  content: "\f13b";
}
.fa-css3:before {
  content: "\f13c";
}
.fa-anchor:before {
  content: "\f13d";
}
.fa-unlock-alt:before {
  content: "\f13e";
}
.fa-bullseye:before {
  content: "\f140";
}
.fa-ellipsis-h:before {
  content: "\f141";
}
.fa-ellipsis-v:before {
  content: "\f142";
}
.fa-rss-square:before {
  content: "\f143";
}
.fa-play-circle:before {
  content: "\f144";
}
.fa-ticket:before {
  content: "\f145";
}
.fa-minus-square:before {
  content: "\f146";
}
.fa-minus-square-o:before {
  content: "\f147";
}
.fa-level-up:before {
  content: "\f148";
}
.fa-level-down:before {
  content: "\f149";
}
.fa-check-square:before {
  content: "\f14a";
}
.fa-pencil-square:before {
  content: "\f14b";
}
.fa-external-link-square:before {
  content: "\f14c";
}
.fa-share-square:before {
  content: "\f14d";
}
.fa-compass:before {
  content: "\f14e";
}
.fa-toggle-down:before,
.fa-caret-square-o-down:before {
  content: "\f150";
}
.fa-toggle-up:before,
.fa-caret-square-o-up:before {
  content: "\f151";
}
.fa-toggle-right:before,
.fa-caret-square-o-right:before {
  content: "\f152";
}
.fa-euro:before,
.fa-eur:before {
  content: "\f153";
}
.fa-gbp:before {
  content: "\f154";
}
.fa-dollar:before,
.fa-usd:before {
  content: "\f155";
}
.fa-rupee:before,
.fa-inr:before {
  content: "\f156";
}
.fa-cny:before,
.fa-rmb:before,
.fa-yen:before,
.fa-jpy:before {
  content: "\f157";
}
.fa-ruble:before,
.fa-rouble:before,
.fa-rub:before {
  content: "\f158";
}
.fa-won:before,
.fa-krw:before {
  content: "\f159";
}
.fa-bitcoin:before,
.fa-btc:before {
  content: "\f15a";
}
.fa-file:before {
  content: "\f15b";
}
.fa-file-text:before {
  content: "\f15c";
}
.fa-sort-alpha-asc:before {
  content: "\f15d";
}
.fa-sort-alpha-desc:before {
  content: "\f15e";
}
.fa-sort-amount-asc:before {
  content: "\f160";
}
.fa-sort-amount-desc:before {
  content: "\f161";
}
.fa-sort-numeric-asc:before {
  content: "\f162";
}
.fa-sort-numeric-desc:before {
  content: "\f163";
}
.fa-thumbs-up:before {
  content: "\f164";
}
.fa-thumbs-down:before {
  content: "\f165";
}
.fa-youtube-square:before {
  content: "\f166";
}
.fa-youtube:before {
  content: "\f167";
}
.fa-xing:before {
  content: "\f168";
}
.fa-xing-square:before {
  content: "\f169";
}
.fa-youtube-play:before {
  content: "\f16a";
}
.fa-dropbox:before {
  content: "\f16b";
}
.fa-stack-overflow:before {
  content: "\f16c";
}
.fa-instagram:before {
  content: "\f16d";
}
.fa-flickr:before {
  content: "\f16e";
}
.fa-adn:before {
  content: "\f170";
}
.fa-bitbucket:before {
  content: "\f171";
}
.fa-bitbucket-square:before {
  content: "\f172";
}
.fa-tumblr:before {
  content: "\f173";
}
.fa-tumblr-square:before {
  content: "\f174";
}
.fa-long-arrow-down:before {
  content: "\f175";
}
.fa-long-arrow-up:before {
  content: "\f176";
}
.fa-long-arrow-left:before {
  content: "\f177";
}
.fa-long-arrow-right:before {
  content: "\f178";
}
.fa-apple:before {
  content: "\f179";
}
.fa-windows:before {
  content: "\f17a";
}
.fa-android:before {
  content: "\f17b";
}
.fa-linux:before {
  content: "\f17c";
}
.fa-dribbble:before {
  content: "\f17d";
}
.fa-skype:before {
  content: "\f17e";
}
.fa-foursquare:before {
  content: "\f180";
}
.fa-trello:before {
  content: "\f181";
}
.fa-female:before {
  content: "\f182";
}
.fa-male:before {
  content: "\f183";
}
.fa-gittip:before {
  content: "\f184";
}
.fa-sun-o:before {
  content: "\f185";
}
.fa-moon-o:before {
  content: "\f186";
}
.fa-archive:before {
  content: "\f187";
}
.fa-bug:before {
  content: "\f188";
}
.fa-vk:before {
  content: "\f189";
}
.fa-weibo:before {
  content: "\f18a";
}
.fa-renren:before {
  content: "\f18b";
}
.fa-pagelines:before {
  content: "\f18c";
}
.fa-stack-exchange:before {
  content: "\f18d";
}
.fa-arrow-circle-o-right:before {
  content: "\f18e";
}
.fa-arrow-circle-o-left:before {
  content: "\f190";
}
.fa-toggle-left:before,
.fa-caret-square-o-left:before {
  content: "\f191";
}
.fa-dot-circle-o:before {
  content: "\f192";
}
.fa-wheelchair:before {
  content: "\f193";
}
.fa-vimeo-square:before {
  content: "\f194";
}
.fa-turkish-lira:before,
.fa-try:before {
  content: "\f195";
}
.fa-plus-square-o:before {
  content: "\f196";
}
.fa-space-shuttle:before {
  content: "\f197";
}
.fa-slack:before {
  content: "\f198";
}
.fa-envelope-square:before {
  content: "\f199";
}
.fa-wordpress:before {
  content: "\f19a";
}
.fa-openid:before {
  content: "\f19b";
}
.fa-institution:before,
.fa-bank:before,
.fa-university:before {
  content: "\f19c";
}
.fa-mortar-board:before,
.fa-graduation-cap:before {
  content: "\f19d";
}
.fa-yahoo:before {
  content: "\f19e";
}
.fa-google:before {
  content: "\f1a0";
}
.fa-reddit:before {
  content: "\f1a1";
}
.fa-reddit-square:before {
  content: "\f1a2";
}
.fa-stumbleupon-circle:before {
  content: "\f1a3";
}
.fa-stumbleupon:before {
  content: "\f1a4";
}
.fa-delicious:before {
  content: "\f1a5";
}
.fa-digg:before {
  content: "\f1a6";
}
.fa-pied-piper:before {
  content: "\f1a7";
}
.fa-pied-piper-alt:before {
  content: "\f1a8";
}
.fa-drupal:before {
  content: "\f1a9";
}
.fa-joomla:before {
  content: "\f1aa";
}
.fa-language:before {
  content: "\f1ab";
}
.fa-fax:before {
  content: "\f1ac";
}
.fa-building:before {
  content: "\f1ad";
}
.fa-child:before {
  content: "\f1ae";
}
.fa-paw:before {
  content: "\f1b0";
}
.fa-spoon:before {
  content: "\f1b1";
}
.fa-cube:before {
  content: "\f1b2";
}
.fa-cubes:before {
  content: "\f1b3";
}
.fa-behance:before {
  content: "\f1b4";
}
.fa-behance-square:before {
  content: "\f1b5";
}
.fa-steam:before {
  content: "\f1b6";
}
.fa-steam-square:before {
  content: "\f1b7";
}
.fa-recycle:before {
  content: "\f1b8";
}
.fa-automobile:before,
.fa-car:before {
  content: "\f1b9";
}
.fa-cab:before,
.fa-taxi:before {
  content: "\f1ba";
}
.fa-tree:before {
  content: "\f1bb";
}
.fa-spotify:before {
  content: "\f1bc";
}
.fa-deviantart:before {
  content: "\f1bd";
}
.fa-soundcloud:before {
  content: "\f1be";
}
.fa-database:before {
  content: "\f1c0";
}
.fa-file-pdf-o:before {
  content: "\f1c1";
}
.fa-file-word-o:before {
  content: "\f1c2";
}
.fa-file-excel-o:before {
  content: "\f1c3";
}
.fa-file-powerpoint-o:before {
  content: "\f1c4";
}
.fa-file-photo-o:before,
.fa-file-picture-o:before,
.fa-file-image-o:before {
  content: "\f1c5";
}
.fa-file-zip-o:before,
.fa-file-archive-o:before {
  content: "\f1c6";
}
.fa-file-sound-o:before,
.fa-file-audio-o:before {
  content: "\f1c7";
}
.fa-file-movie-o:before,
.fa-file-video-o:before {
  content: "\f1c8";
}
.fa-file-code-o:before {
  content: "\f1c9";
}
.fa-vine:before {
  content: "\f1ca";
}
.fa-codepen:before {
  content: "\f1cb";
}
.fa-jsfiddle:before {
  content: "\f1cc";
}
.fa-life-bouy:before,
.fa-life-buoy:before,
.fa-life-saver:before,
.fa-support:before,
.fa-life-ring:before {
  content: "\f1cd";
}
.fa-circle-o-notch:before {
  content: "\f1ce";
}
.fa-ra:before,
.fa-rebel:before {
  content: "\f1d0";
}
.fa-ge:before,
.fa-empire:before {
  content: "\f1d1";
}
.fa-git-square:before {
  content: "\f1d2";
}
.fa-git:before {
  content: "\f1d3";
}
.fa-hacker-news:before {
  content: "\f1d4";
}
.fa-tencent-weibo:before {
  content: "\f1d5";
}
.fa-qq:before {
  content: "\f1d6";
}
.fa-wechat:before,
.fa-weixin:before {
  content: "\f1d7";
}
.fa-send:before,
.fa-paper-plane:before {
  content: "\f1d8";
}
.fa-send-o:before,
.fa-paper-plane-o:before {
  content: "\f1d9";
}
.fa-history:before {
  content: "\f1da";
}
.fa-circle-thin:before {
  content: "\f1db";
}
.fa-header:before {
  content: "\f1dc";
}
.fa-paragraph:before {
  content: "\f1dd";
}
.fa-sliders:before {
  content: "\f1de";
}
.fa-share-alt:before {
  content: "\f1e0";
}
.fa-share-alt-square:before {
  content: "\f1e1";
}
.fa-bomb:before {
  content: "\f1e2";
}
.fa-soccer-ball-o:before,
.fa-futbol-o:before {
  content: "\f1e3";
}
.fa-tty:before {
  content: "\f1e4";
}
.fa-binoculars:before {
  content: "\f1e5";
}
.fa-plug:before {
  content: "\f1e6";
}
.fa-slideshare:before {
  content: "\f1e7";
}
.fa-twitch:before {
  content: "\f1e8";
}
.fa-yelp:before {
  content: "\f1e9";
}
.fa-newspaper-o:before {
  content: "\f1ea";
}
.fa-wifi:before {
  content: "\f1eb";
}
.fa-calculator:before {
  content: "\f1ec";
}
.fa-paypal:before {
  content: "\f1ed";
}
.fa-google-wallet:before {
  content: "\f1ee";
}
.fa-cc-visa:before {
  content: "\f1f0";
}
.fa-cc-mastercard:before {
  content: "\f1f1";
}
.fa-cc-discover:before {
  content: "\f1f2";
}
.fa-cc-amex:before {
  content: "\f1f3";
}
.fa-cc-paypal:before {
  content: "\f1f4";
}
.fa-cc-stripe:before {
  content: "\f1f5";
}
.fa-bell-slash:before {
  content: "\f1f6";
}
.fa-bell-slash-o:before {
  content: "\f1f7";
}
.fa-trash:before {
  content: "\f1f8";
}
.fa-copyright:before {
  content: "\f1f9";
}
.fa-at:before {
  content: "\f1fa";
}
.fa-eyedropper:before {
  content: "\f1fb";
}
.fa-paint-brush:before {
  content: "\f1fc";
}
.fa-birthday-cake:before {
  content: "\f1fd";
}
.fa-area-chart:before {
  content: "\f1fe";
}
.fa-pie-chart:before {
  content: "\f200";
}
.fa-line-chart:before {
  content: "\f201";
}
.fa-lastfm:before {
  content: "\f202";
}
.fa-lastfm-square:before {
  content: "\f203";
}
.fa-toggle-off:before {
  content: "\f204";
}
.fa-toggle-on:before {
  content: "\f205";
}
.fa-bicycle:before {
  content: "\f206";
}
.fa-bus:before {
  content: "\f207";
}
.fa-ioxhost:before {
  content: "\f208";
}
.fa-angellist:before {
  content: "\f209";
}
.fa-cc:before {
  content: "\f20a";
}
.fa-shekel:before,
.fa-sheqel:before,
.fa-ils:before {
  content: "\f20b";
}
.fa-meanpath:before {
  content: "\f20c";
}
/*!
*
* IPython base
*
*/
.modal.fade .modal-dialog {
  -webkit-transform: translate(0, 0);
  -ms-transform: translate(0, 0);
  -o-transform: translate(0, 0);
  transform: translate(0, 0);
}
code {
  color: #000;
}
pre {
  font-size: inherit;
  line-height: inherit;
}
label {
  font-weight: normal;
}
/* Make the page background atleast 100% the height of the view port */
/* Make the page itself atleast 70% the height of the view port */
.border-box-sizing {
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}
.corner-all {
  border-radius: 2px;
}
.no-padding {
  padding: 0px;
}
/* Flexible box model classes */
/* Taken from Alex Russell http://infrequently.org/2009/08/css-3-progress/ */
/* This file is a compatability layer.  It allows the usage of flexible box 
model layouts accross multiple browsers, including older browsers.  The newest,
universal implementation of the flexible box model is used when available (see
`Modern browsers` comments below).  Browsers that are known to implement this 
new spec completely include:

    Firefox 28.0+
    Chrome 29.0+
    Internet Explorer 11+ 
    Opera 17.0+

Browsers not listed, including Safari, are supported via the styling under the
`Old browsers` comments below.
*/
.hbox {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
}
.hbox > * {
  /* Old browsers */
  -webkit-box-flex: 0;
  -moz-box-flex: 0;
  box-flex: 0;
  /* Modern browsers */
  flex: none;
}
.vbox {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
}
.vbox > * {
  /* Old browsers */
  -webkit-box-flex: 0;
  -moz-box-flex: 0;
  box-flex: 0;
  /* Modern browsers */
  flex: none;
}
.hbox.reverse,
.vbox.reverse,
.reverse {
  /* Old browsers */
  -webkit-box-direction: reverse;
  -moz-box-direction: reverse;
  box-direction: reverse;
  /* Modern browsers */
  flex-direction: row-reverse;
}
.hbox.box-flex0,
.vbox.box-flex0,
.box-flex0 {
  /* Old browsers */
  -webkit-box-flex: 0;
  -moz-box-flex: 0;
  box-flex: 0;
  /* Modern browsers */
  flex: none;
  width: auto;
}
.hbox.box-flex1,
.vbox.box-flex1,
.box-flex1 {
  /* Old browsers */
  -webkit-box-flex: 1;
  -moz-box-flex: 1;
  box-flex: 1;
  /* Modern browsers */
  flex: 1;
}
.hbox.box-flex,
.vbox.box-flex,
.box-flex {
  /* Old browsers */
  /* Old browsers */
  -webkit-box-flex: 1;
  -moz-box-flex: 1;
  box-flex: 1;
  /* Modern browsers */
  flex: 1;
}
.hbox.box-flex2,
.vbox.box-flex2,
.box-flex2 {
  /* Old browsers */
  -webkit-box-flex: 2;
  -moz-box-flex: 2;
  box-flex: 2;
  /* Modern browsers */
  flex: 2;
}
.box-group1 {
  /*  Deprecated */
  -webkit-box-flex-group: 1;
  -moz-box-flex-group: 1;
  box-flex-group: 1;
}
.box-group2 {
  /* Deprecated */
  -webkit-box-flex-group: 2;
  -moz-box-flex-group: 2;
  box-flex-group: 2;
}
.hbox.start,
.vbox.start,
.start {
  /* Old browsers */
  -webkit-box-pack: start;
  -moz-box-pack: start;
  box-pack: start;
  /* Modern browsers */
  justify-content: flex-start;
}
.hbox.end,
.vbox.end,
.end {
  /* Old browsers */
  -webkit-box-pack: end;
  -moz-box-pack: end;
  box-pack: end;
  /* Modern browsers */
  justify-content: flex-end;
}
.hbox.center,
.vbox.center,
.center {
  /* Old browsers */
  -webkit-box-pack: center;
  -moz-box-pack: center;
  box-pack: center;
  /* Modern browsers */
  justify-content: center;
}
.hbox.baseline,
.vbox.baseline,
.baseline {
  /* Old browsers */
  -webkit-box-pack: baseline;
  -moz-box-pack: baseline;
  box-pack: baseline;
  /* Modern browsers */
  justify-content: baseline;
}
.hbox.stretch,
.vbox.stretch,
.stretch {
  /* Old browsers */
  -webkit-box-pack: stretch;
  -moz-box-pack: stretch;
  box-pack: stretch;
  /* Modern browsers */
  justify-content: stretch;
}
.hbox.align-start,
.vbox.align-start,
.align-start {
  /* Old browsers */
  -webkit-box-align: start;
  -moz-box-align: start;
  box-align: start;
  /* Modern browsers */
  align-items: flex-start;
}
.hbox.align-end,
.vbox.align-end,
.align-end {
  /* Old browsers */
  -webkit-box-align: end;
  -moz-box-align: end;
  box-align: end;
  /* Modern browsers */
  align-items: flex-end;
}
.hbox.align-center,
.vbox.align-center,
.align-center {
  /* Old browsers */
  -webkit-box-align: center;
  -moz-box-align: center;
  box-align: center;
  /* Modern browsers */
  align-items: center;
}
.hbox.align-baseline,
.vbox.align-baseline,
.align-baseline {
  /* Old browsers */
  -webkit-box-align: baseline;
  -moz-box-align: baseline;
  box-align: baseline;
  /* Modern browsers */
  align-items: baseline;
}
.hbox.align-stretch,
.vbox.align-stretch,
.align-stretch {
  /* Old browsers */
  -webkit-box-align: stretch;
  -moz-box-align: stretch;
  box-align: stretch;
  /* Modern browsers */
  align-items: stretch;
}
div.error {
  margin: 2em;
  text-align: center;
}
div.error > h1 {
  font-size: 500%;
  line-height: normal;
}
div.error > p {
  font-size: 200%;
  line-height: normal;
}
div.traceback-wrapper {
  text-align: left;
  max-width: 800px;
  margin: auto;
}
/**
 * Primary styles
 *
 * Author: Jupyter Development Team
 */
body {
  background-color: #fff;
  /* This makes sure that the body covers the entire window and needs to
       be in a different element than the display: box in wrapper below */
  position: absolute;
  left: 0px;
  right: 0px;
  top: 0px;
  bottom: 0px;
  overflow: visible;
}
body > #header {
  /* Initially hidden to prevent FLOUC */
  display: none;
  background-color: #fff;
  /* Display over codemirror */
  position: relative;
  z-index: 100;
}
body > #header #header-container {
  padding-bottom: 5px;
  padding-top: 5px;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}
body > #header .header-bar {
  width: 100%;
  height: 1px;
  background: #e7e7e7;
  margin-bottom: -1px;
}
@media print {
  body > #header {
    display: none !important;
  }
}
#header-spacer {
  width: 100%;
  visibility: hidden;
}
@media print {
  #header-spacer {
    display: none;
  }
}
#ipython_notebook {
  padding-left: 0px;
  padding-top: 1px;
  padding-bottom: 1px;
}
@media (max-width: 991px) {
  #ipython_notebook {
    margin-left: 10px;
  }
}
[dir="rtl"] #ipython_notebook {
  float: right !important;
}
#noscript {
  width: auto;
  padding-top: 16px;
  padding-bottom: 16px;
  text-align: center;
  font-size: 22px;
  color: red;
  font-weight: bold;
}
#ipython_notebook img {
  height: 28px;
}
#site {
  width: 100%;
  display: none;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  overflow: auto;
}
@media print {
  #site {
    height: auto !important;
  }
}
/* Smaller buttons */
.ui-button .ui-button-text {
  padding: 0.2em 0.8em;
  font-size: 77%;
}
input.ui-button {
  padding: 0.3em 0.9em;
}
span#login_widget {
  float: right;
}
span#login_widget > .button,
#logout {
  color: #333;
  background-color: #fff;
  border-color: #ccc;
}
span#login_widget > .button:focus,
#logout:focus,
span#login_widget > .button.focus,
#logout.focus {
  color: #333;
  background-color: #e6e6e6;
  border-color: #8c8c8c;
}
span#login_widget > .button:hover,
#logout:hover {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
span#login_widget > .button:active,
#logout:active,
span#login_widget > .button.active,
#logout.active,
.open > .dropdown-togglespan#login_widget > .button,
.open > .dropdown-toggle#logout {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
span#login_widget > .button:active:hover,
#logout:active:hover,
span#login_widget > .button.active:hover,
#logout.active:hover,
.open > .dropdown-togglespan#login_widget > .button:hover,
.open > .dropdown-toggle#logout:hover,
span#login_widget > .button:active:focus,
#logout:active:focus,
span#login_widget > .button.active:focus,
#logout.active:focus,
.open > .dropdown-togglespan#login_widget > .button:focus,
.open > .dropdown-toggle#logout:focus,
span#login_widget > .button:active.focus,
#logout:active.focus,
span#login_widget > .button.active.focus,
#logout.active.focus,
.open > .dropdown-togglespan#login_widget > .button.focus,
.open > .dropdown-toggle#logout.focus {
  color: #333;
  background-color: #d4d4d4;
  border-color: #8c8c8c;
}
span#login_widget > .button:active,
#logout:active,
span#login_widget > .button.active,
#logout.active,
.open > .dropdown-togglespan#login_widget > .button,
.open > .dropdown-toggle#logout {
  background-image: none;
}
span#login_widget > .button.disabled:hover,
#logout.disabled:hover,
span#login_widget > .button[disabled]:hover,
#logout[disabled]:hover,
fieldset[disabled] span#login_widget > .button:hover,
fieldset[disabled] #logout:hover,
span#login_widget > .button.disabled:focus,
#logout.disabled:focus,
span#login_widget > .button[disabled]:focus,
#logout[disabled]:focus,
fieldset[disabled] span#login_widget > .button:focus,
fieldset[disabled] #logout:focus,
span#login_widget > .button.disabled.focus,
#logout.disabled.focus,
span#login_widget > .button[disabled].focus,
#logout[disabled].focus,
fieldset[disabled] span#login_widget > .button.focus,
fieldset[disabled] #logout.focus {
  background-color: #fff;
  border-color: #ccc;
}
span#login_widget > .button .badge,
#logout .badge {
  color: #fff;
  background-color: #333;
}
.nav-header {
  text-transform: none;
}
#header > span {
  margin-top: 10px;
}
.modal_stretch .modal-dialog {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
  min-height: 80vh;
}
.modal_stretch .modal-dialog .modal-body {
  max-height: calc(100vh - 200px);
  overflow: auto;
  flex: 1;
}
@media (min-width: 768px) {
  .modal .modal-dialog {
    width: 700px;
  }
}
@media (min-width: 768px) {
  select.form-control {
    margin-left: 12px;
    margin-right: 12px;
  }
}
/*!
*
* IPython auth
*
*/
.center-nav {
  display: inline-block;
  margin-bottom: -4px;
}
/*!
*
* IPython tree view
*
*/
/* We need an invisible input field on top of the sentense*/
/* "Drag file onto the list ..." */
.alternate_upload {
  background-color: none;
  display: inline;
}
.alternate_upload.form {
  padding: 0;
  margin: 0;
}
.alternate_upload input.fileinput {
  text-align: center;
  vertical-align: middle;
  display: inline;
  opacity: 0;
  z-index: 2;
  width: 12ex;
  margin-right: -12ex;
}
.alternate_upload .btn-upload {
  height: 22px;
}
/**
 * Primary styles
 *
 * Author: Jupyter Development Team
 */
[dir="rtl"] #tabs li {
  float: right;
}
ul#tabs {
  margin-bottom: 4px;
}
[dir="rtl"] ul#tabs {
  margin-right: 0px;
}
ul#tabs a {
  padding-top: 6px;
  padding-bottom: 4px;
}
ul.breadcrumb a:focus,
ul.breadcrumb a:hover {
  text-decoration: none;
}
ul.breadcrumb i.icon-home {
  font-size: 16px;
  margin-right: 4px;
}
ul.breadcrumb span {
  color: #5e5e5e;
}
.list_toolbar {
  padding: 4px 0 4px 0;
  vertical-align: middle;
}
.list_toolbar .tree-buttons {
  padding-top: 1px;
}
[dir="rtl"] .list_toolbar .tree-buttons {
  float: left !important;
}
[dir="rtl"] .list_toolbar .pull-right {
  padding-top: 1px;
  float: left !important;
}
[dir="rtl"] .list_toolbar .pull-left {
  float: right !important;
}
.dynamic-buttons {
  padding-top: 3px;
  display: inline-block;
}
.list_toolbar [class*="span"] {
  min-height: 24px;
}
.list_header {
  font-weight: bold;
  background-color: #EEE;
}
.list_placeholder {
  font-weight: bold;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 7px;
  padding-right: 7px;
}
.list_container {
  margin-top: 4px;
  margin-bottom: 20px;
  border: 1px solid #ddd;
  border-radius: 2px;
}
.list_container > div {
  border-bottom: 1px solid #ddd;
}
.list_container > div:hover .list-item {
  background-color: red;
}
.list_container > div:last-child {
  border: none;
}
.list_item:hover .list_item {
  background-color: #ddd;
}
.list_item a {
  text-decoration: none;
}
.list_item:hover {
  background-color: #fafafa;
}
.list_header > div,
.list_item > div {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 7px;
  padding-right: 7px;
  line-height: 22px;
}
.list_header > div input,
.list_item > div input {
  margin-right: 7px;
  margin-left: 14px;
  vertical-align: baseline;
  line-height: 22px;
  position: relative;
  top: -1px;
}
.list_header > div .item_link,
.list_item > div .item_link {
  margin-left: -1px;
  vertical-align: baseline;
  line-height: 22px;
}
.new-file input[type=checkbox] {
  visibility: hidden;
}
.item_name {
  line-height: 22px;
  height: 24px;
}
.item_icon {
  font-size: 14px;
  color: #5e5e5e;
  margin-right: 7px;
  margin-left: 7px;
  line-height: 22px;
  vertical-align: baseline;
}
.item_buttons {
  line-height: 1em;
  margin-left: -5px;
}
.item_buttons .btn,
.item_buttons .btn-group,
.item_buttons .input-group {
  float: left;
}
.item_buttons > .btn,
.item_buttons > .btn-group,
.item_buttons > .input-group {
  margin-left: 5px;
}
.item_buttons .btn {
  min-width: 13ex;
}
.item_buttons .running-indicator {
  padding-top: 4px;
  color: #5cb85c;
}
.item_buttons .kernel-name {
  padding-top: 4px;
  color: #5bc0de;
  margin-right: 7px;
  float: left;
}
.toolbar_info {
  height: 24px;
  line-height: 24px;
}
.list_item input:not([type=checkbox]) {
  padding-top: 3px;
  padding-bottom: 3px;
  height: 22px;
  line-height: 14px;
  margin: 0px;
}
.highlight_text {
  color: blue;
}
#project_name {
  display: inline-block;
  padding-left: 7px;
  margin-left: -2px;
}
#project_name > .breadcrumb {
  padding: 0px;
  margin-bottom: 0px;
  background-color: transparent;
  font-weight: bold;
}
#tree-selector {
  padding-right: 0px;
}
[dir="rtl"] #tree-selector a {
  float: right;
}
#button-select-all {
  min-width: 50px;
}
#select-all {
  margin-left: 7px;
  margin-right: 2px;
}
.menu_icon {
  margin-right: 2px;
}
.tab-content .row {
  margin-left: 0px;
  margin-right: 0px;
}
.folder_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f114";
}
.folder_icon:before.pull-left {
  margin-right: .3em;
}
.folder_icon:before.pull-right {
  margin-left: .3em;
}
.notebook_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f02d";
  position: relative;
  top: -1px;
}
.notebook_icon:before.pull-left {
  margin-right: .3em;
}
.notebook_icon:before.pull-right {
  margin-left: .3em;
}
.running_notebook_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f02d";
  position: relative;
  top: -1px;
  color: #5cb85c;
}
.running_notebook_icon:before.pull-left {
  margin-right: .3em;
}
.running_notebook_icon:before.pull-right {
  margin-left: .3em;
}
.file_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f016";
  position: relative;
  top: -2px;
}
.file_icon:before.pull-left {
  margin-right: .3em;
}
.file_icon:before.pull-right {
  margin-left: .3em;
}
#notebook_toolbar .pull-right {
  padding-top: 0px;
  margin-right: -1px;
}
ul#new-menu {
  left: auto;
  right: 0;
}
[dir="rtl"] #new-menu {
  text-align: right;
}
.kernel-menu-icon {
  padding-right: 12px;
  width: 24px;
  content: "\f096";
}
.kernel-menu-icon:before {
  content: "\f096";
}
.kernel-menu-icon-current:before {
  content: "\f00c";
}
#tab_content {
  padding-top: 20px;
}
#running .panel-group .panel {
  margin-top: 3px;
  margin-bottom: 1em;
}
#running .panel-group .panel .panel-heading {
  background-color: #EEE;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 7px;
  padding-right: 7px;
  line-height: 22px;
}
#running .panel-group .panel .panel-heading a:focus,
#running .panel-group .panel .panel-heading a:hover {
  text-decoration: none;
}
#running .panel-group .panel .panel-body {
  padding: 0px;
}
#running .panel-group .panel .panel-body .list_container {
  margin-top: 0px;
  margin-bottom: 0px;
  border: 0px;
  border-radius: 0px;
}
#running .panel-group .panel .panel-body .list_container .list_item {
  border-bottom: 1px solid #ddd;
}
#running .panel-group .panel .panel-body .list_container .list_item:last-child {
  border-bottom: 0px;
}
[dir="rtl"] #running .col-sm-8 {
  float: right !important;
}
.delete-button {
  display: none;
}
.duplicate-button {
  display: none;
}
.rename-button {
  display: none;
}
.shutdown-button {
  display: none;
}
.dynamic-instructions {
  display: inline-block;
  padding-top: 4px;
}
/*!
*
* IPython text editor webapp
*
*/
.selected-keymap i.fa {
  padding: 0px 5px;
}
.selected-keymap i.fa:before {
  content: "\f00c";
}
#mode-menu {
  overflow: auto;
  max-height: 20em;
}
.edit_app #header {
  -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
}
.edit_app #menubar .navbar {
  /* Use a negative 1 bottom margin, so the border overlaps the border of the
    header */
  margin-bottom: -1px;
}
.dirty-indicator {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  width: 20px;
}
.dirty-indicator.pull-left {
  margin-right: .3em;
}
.dirty-indicator.pull-right {
  margin-left: .3em;
}
.dirty-indicator-dirty {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  width: 20px;
}
.dirty-indicator-dirty.pull-left {
  margin-right: .3em;
}
.dirty-indicator-dirty.pull-right {
  margin-left: .3em;
}
.dirty-indicator-clean {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  width: 20px;
}
.dirty-indicator-clean.pull-left {
  margin-right: .3em;
}
.dirty-indicator-clean.pull-right {
  margin-left: .3em;
}
.dirty-indicator-clean:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f00c";
}
.dirty-indicator-clean:before.pull-left {
  margin-right: .3em;
}
.dirty-indicator-clean:before.pull-right {
  margin-left: .3em;
}
#filename {
  font-size: 16pt;
  display: table;
  padding: 0px 5px;
}
#current-mode {
  padding-left: 5px;
  padding-right: 5px;
}
#texteditor-backdrop {
  padding-top: 20px;
  padding-bottom: 20px;
}
@media not print {
  #texteditor-backdrop {
    background-color: #EEE;
  }
}
@media print {
  #texteditor-backdrop #texteditor-container .CodeMirror-gutter,
  #texteditor-backdrop #texteditor-container .CodeMirror-gutters {
    background-color: #fff;
  }
}
@media not print {
  #texteditor-backdrop #texteditor-container .CodeMirror-gutter,
  #texteditor-backdrop #texteditor-container .CodeMirror-gutters {
    background-color: #fff;
  }
}
@media not print {
  #texteditor-backdrop #texteditor-container {
    padding: 0px;
    background-color: #fff;
    -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
    box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  }
}
/*!
*
* IPython notebook
*
*/
/* CSS font colors for translated ANSI colors. */
.ansibold {
  font-weight: bold;
}
/* use dark versions for foreground, to improve visibility */
.ansiblack {
  color: black;
}
.ansired {
  color: darkred;
}
.ansigreen {
  color: darkgreen;
}
.ansiyellow {
  color: #c4a000;
}
.ansiblue {
  color: darkblue;
}
.ansipurple {
  color: darkviolet;
}
.ansicyan {
  color: steelblue;
}
.ansigray {
  color: gray;
}
/* and light for background, for the same reason */
.ansibgblack {
  background-color: black;
}
.ansibgred {
  background-color: red;
}
.ansibggreen {
  background-color: green;
}
.ansibgyellow {
  background-color: yellow;
}
.ansibgblue {
  background-color: blue;
}
.ansibgpurple {
  background-color: magenta;
}
.ansibgcyan {
  background-color: cyan;
}
.ansibggray {
  background-color: gray;
}
div.cell {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
  border-radius: 2px;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  border-width: 1px;
  border-style: solid;
  border-color: transparent;
  width: 100%;
  padding: 5px;
  /* This acts as a spacer between cells, that is outside the border */
  margin: 0px;
  outline: none;
  border-left-width: 1px;
  padding-left: 5px;
  background: linear-gradient(to right, transparent -40px, transparent 1px, transparent 1px, transparent 100%);
}
div.cell.jupyter-soft-selected {
  border-left-color: #90CAF9;
  border-left-color: #E3F2FD;
  border-left-width: 1px;
  padding-left: 5px;
  border-right-color: #E3F2FD;
  border-right-width: 1px;
  background: #E3F2FD;
}
@media print {
  div.cell.jupyter-soft-selected {
    border-color: transparent;
  }
}
div.cell.selected {
  border-color: #ababab;
  border-left-width: 0px;
  padding-left: 6px;
  background: linear-gradient(to right, #42A5F5 -40px, #42A5F5 5px, transparent 5px, transparent 100%);
}
@media print {
  div.cell.selected {
    border-color: transparent;
  }
}
div.cell.selected.jupyter-soft-selected {
  border-left-width: 0;
  padding-left: 6px;
  background: linear-gradient(to right, #42A5F5 -40px, #42A5F5 7px, #E3F2FD 7px, #E3F2FD 100%);
}
.edit_mode div.cell.selected {
  border-color: #66BB6A;
  border-left-width: 0px;
  padding-left: 6px;
  background: linear-gradient(to right, #66BB6A -40px, #66BB6A 5px, transparent 5px, transparent 100%);
}
@media print {
  .edit_mode div.cell.selected {
    border-color: transparent;
  }
}
.prompt {
  /* This needs to be wide enough for 3 digit prompt numbers: In[100]: */
  min-width: 14ex;
  /* This padding is tuned to match the padding on the CodeMirror editor. */
  padding: 0.4em;
  margin: 0px;
  font-family: monospace;
  text-align: right;
  /* This has to match that of the the CodeMirror class line-height below */
  line-height: 1.21429em;
  /* Don't highlight prompt number selection */
  -webkit-touch-callout: none;
  -webkit-user-select: none;
  -khtml-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
  /* Use default cursor */
  cursor: default;
}
@media (max-width: 540px) {
  .prompt {
    text-align: left;
  }
}
div.inner_cell {
  min-width: 0;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
  /* Old browsers */
  -webkit-box-flex: 1;
  -moz-box-flex: 1;
  box-flex: 1;
  /* Modern browsers */
  flex: 1;
}
/* input_area and input_prompt must match in top border and margin for alignment */
div.input_area {
  border: 1px solid #cfcfcf;
  border-radius: 2px;
  background: #f7f7f7;
  line-height: 1.21429em;
}
/* This is needed so that empty prompt areas can collapse to zero height when there
   is no content in the output_subarea and the prompt. The main purpose of this is
   to make sure that empty JavaScript output_subareas have no height. */
div.prompt:empty {
  padding-top: 0;
  padding-bottom: 0;
}
div.unrecognized_cell {
  padding: 5px 5px 5px 0px;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
}
div.unrecognized_cell .inner_cell {
  border-radius: 2px;
  padding: 5px;
  font-weight: bold;
  color: red;
  border: 1px solid #cfcfcf;
  background: #eaeaea;
}
div.unrecognized_cell .inner_cell a {
  color: inherit;
  text-decoration: none;
}
div.unrecognized_cell .inner_cell a:hover {
  color: inherit;
  text-decoration: none;
}
@media (max-width: 540px) {
  div.unrecognized_cell > div.prompt {
    display: none;
  }
}
div.code_cell {
  /* avoid page breaking on code cells when printing */
}
@media print {
  div.code_cell {
    page-break-inside: avoid;
  }
}
/* any special styling for code cells that are currently running goes here */
div.input {
  page-break-inside: avoid;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
}
@media (max-width: 540px) {
  div.input {
    /* Old browsers */
    display: -webkit-box;
    -webkit-box-orient: vertical;
    -webkit-box-align: stretch;
    display: -moz-box;
    -moz-box-orient: vertical;
    -moz-box-align: stretch;
    display: box;
    box-orient: vertical;
    box-align: stretch;
    /* Modern browsers */
    display: flex;
    flex-direction: column;
    align-items: stretch;
  }
}
/* input_area and input_prompt must match in top border and margin for alignment */
div.input_prompt {
  color: #303F9F;
  border-top: 1px solid transparent;
}
div.input_area > div.highlight {
  margin: 0.4em;
  border: none;
  padding: 0px;
  background-color: transparent;
}
div.input_area > div.highlight > pre {
  margin: 0px;
  border: none;
  padding: 0px;
  background-color: transparent;
}
/* The following gets added to the <head> if it is detected that the user has a
 * monospace font with inconsistent normal/bold/italic height.  See
 * notebookmain.js.  Such fonts will have keywords vertically offset with
 * respect to the rest of the text.  The user should select a better font.
 * See: https://github.com/ipython/ipython/issues/1503
 *
 * .CodeMirror span {
 *      vertical-align: bottom;
 * }
 */
.CodeMirror {
  line-height: 1.21429em;
  /* Changed from 1em to our global default */
  font-size: 14px;
  height: auto;
  /* Changed to auto to autogrow */
  background: none;
  /* Changed from white to allow our bg to show through */
}
.CodeMirror-scroll {
  /*  The CodeMirror docs are a bit fuzzy on if overflow-y should be hidden or visible.*/
  /*  We have found that if it is visible, vertical scrollbars appear with font size changes.*/
  overflow-y: hidden;
  overflow-x: auto;
}
.CodeMirror-lines {
  /* In CM2, this used to be 0.4em, but in CM3 it went to 4px. We need the em value because */
  /* we have set a different line-height and want this to scale with that. */
  padding: 0.4em;
}
.CodeMirror-linenumber {
  padding: 0 8px 0 4px;
}
.CodeMirror-gutters {
  border-bottom-left-radius: 2px;
  border-top-left-radius: 2px;
}
.CodeMirror pre {
  /* In CM3 this went to 4px from 0 in CM2. We need the 0 value because of how we size */
  /* .CodeMirror-lines */
  padding: 0;
  border: 0;
  border-radius: 0;
}
/*

Original style from softwaremaniacs.org (c) Ivan Sagalaev <Maniac@SoftwareManiacs.Org>
Adapted from GitHub theme

*/
.highlight-base {
  color: #000;
}
.highlight-variable {
  color: #000;
}
.highlight-variable-2 {
  color: #1a1a1a;
}
.highlight-variable-3 {
  color: #333333;
}
.highlight-string {
  color: #BA2121;
}
.highlight-comment {
  color: #408080;
  font-style: italic;
}
.highlight-number {
  color: #080;
}
.highlight-atom {
  color: #88F;
}
.highlight-keyword {
  color: #008000;
  font-weight: bold;
}
.highlight-builtin {
  color: #008000;
}
.highlight-error {
  color: #f00;
}
.highlight-operator {
  color: #AA22FF;
  font-weight: bold;
}
.highlight-meta {
  color: #AA22FF;
}
/* previously not defined, copying from default codemirror */
.highlight-def {
  color: #00f;
}
.highlight-string-2 {
  color: #f50;
}
.highlight-qualifier {
  color: #555;
}
.highlight-bracket {
  color: #997;
}
.highlight-tag {
  color: #170;
}
.highlight-attribute {
  color: #00c;
}
.highlight-header {
  color: blue;
}
.highlight-quote {
  color: #090;
}
.highlight-link {
  color: #00c;
}
/* apply the same style to codemirror */
.cm-s-ipython span.cm-keyword {
  color: #008000;
  font-weight: bold;
}
.cm-s-ipython span.cm-atom {
  color: #88F;
}
.cm-s-ipython span.cm-number {
  color: #080;
}
.cm-s-ipython span.cm-def {
  color: #00f;
}
.cm-s-ipython span.cm-variable {
  color: #000;
}
.cm-s-ipython span.cm-operator {
  color: #AA22FF;
  font-weight: bold;
}
.cm-s-ipython span.cm-variable-2 {
  color: #1a1a1a;
}
.cm-s-ipython span.cm-variable-3 {
  color: #333333;
}
.cm-s-ipython span.cm-comment {
  color: #408080;
  font-style: italic;
}
.cm-s-ipython span.cm-string {
  color: #BA2121;
}
.cm-s-ipython span.cm-string-2 {
  color: #f50;
}
.cm-s-ipython span.cm-meta {
  color: #AA22FF;
}
.cm-s-ipython span.cm-qualifier {
  color: #555;
}
.cm-s-ipython span.cm-builtin {
  color: #008000;
}
.cm-s-ipython span.cm-bracket {
  color: #997;
}
.cm-s-ipython span.cm-tag {
  color: #170;
}
.cm-s-ipython span.cm-attribute {
  color: #00c;
}
.cm-s-ipython span.cm-header {
  color: blue;
}
.cm-s-ipython span.cm-quote {
  color: #090;
}
.cm-s-ipython span.cm-link {
  color: #00c;
}
.cm-s-ipython span.cm-error {
  color: #f00;
}
.cm-s-ipython span.cm-tab {
  background: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADAAAAAMCAYAAAAkuj5RAAAAAXNSR0IArs4c6QAAAGFJREFUSMft1LsRQFAQheHPowAKoACx3IgEKtaEHujDjORSgWTH/ZOdnZOcM/sgk/kFFWY0qV8foQwS4MKBCS3qR6ixBJvElOobYAtivseIE120FaowJPN75GMu8j/LfMwNjh4HUpwg4LUAAAAASUVORK5CYII=);
  background-position: right;
  background-repeat: no-repeat;
}
div.output_wrapper {
  /* this position must be relative to enable descendents to be absolute within it */
  position: relative;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
  z-index: 1;
}
/* class for the output area when it should be height-limited */
div.output_scroll {
  /* ideally, this would be max-height, but FF barfs all over that */
  height: 24em;
  /* FF needs this *and the wrapper* to specify full width, or it will shrinkwrap */
  width: 100%;
  overflow: auto;
  border-radius: 2px;
  -webkit-box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.8);
  box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.8);
  display: block;
}
/* output div while it is collapsed */
div.output_collapsed {
  margin: 0px;
  padding: 0px;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
}
div.out_prompt_overlay {
  height: 100%;
  padding: 0px 0.4em;
  position: absolute;
  border-radius: 2px;
}
div.out_prompt_overlay:hover {
  /* use inner shadow to get border that is computed the same on WebKit/FF */
  -webkit-box-shadow: inset 0 0 1px #000;
  box-shadow: inset 0 0 1px #000;
  background: rgba(240, 240, 240, 0.5);
}
div.output_prompt {
  color: #D84315;
}
/* This class is the outer container of all output sections. */
div.output_area {
  padding: 0px;
  page-break-inside: avoid;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
}
div.output_area .MathJax_Display {
  text-align: left !important;
}
div.output_area .rendered_html table {
  margin-left: 0;
  margin-right: 0;
}
div.output_area .rendered_html img {
  margin-left: 0;
  margin-right: 0;
}
div.output_area img,
div.output_area svg {
  max-width: 100%;
  height: auto;
}
div.output_area img.unconfined,
div.output_area svg.unconfined {
  max-width: none;
}
/* This is needed to protect the pre formating from global settings such
   as that of bootstrap */
.output {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: vertical;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: vertical;
  -moz-box-align: stretch;
  display: box;
  box-orient: vertical;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: column;
  align-items: stretch;
}
@media (max-width: 540px) {
  div.output_area {
    /* Old browsers */
    display: -webkit-box;
    -webkit-box-orient: vertical;
    -webkit-box-align: stretch;
    display: -moz-box;
    -moz-box-orient: vertical;
    -moz-box-align: stretch;
    display: box;
    box-orient: vertical;
    box-align: stretch;
    /* Modern browsers */
    display: flex;
    flex-direction: column;
    align-items: stretch;
  }
}
div.output_area pre {
  margin: 0;
  padding: 0;
  border: 0;
  vertical-align: baseline;
  color: black;
  background-color: transparent;
  border-radius: 0;
}
/* This class is for the output subarea inside the output_area and after
   the prompt div. */
div.output_subarea {
  overflow-x: auto;
  padding: 0.4em;
  /* Old browsers */
  -webkit-box-flex: 1;
  -moz-box-flex: 1;
  box-flex: 1;
  /* Modern browsers */
  flex: 1;
  max-width: calc(100% - 14ex);
}
div.output_scroll div.output_subarea {
  overflow-x: visible;
}
/* The rest of the output_* classes are for special styling of the different
   output types */
/* all text output has this class: */
div.output_text {
  text-align: left;
  color: #000;
  /* This has to match that of the the CodeMirror class line-height below */
  line-height: 1.21429em;
}
/* stdout/stderr are 'text' as well as 'stream', but execute_result/error are *not* streams */
div.output_stderr {
  background: #fdd;
  /* very light red background for stderr */
}
div.output_latex {
  text-align: left;
}
/* Empty output_javascript divs should have no height */
div.output_javascript:empty {
  padding: 0;
}
.js-error {
  color: darkred;
}
/* raw_input styles */
div.raw_input_container {
  line-height: 1.21429em;
  padding-top: 5px;
}
pre.raw_input_prompt {
  /* nothing needed here. */
}
input.raw_input {
  font-family: monospace;
  font-size: inherit;
  color: inherit;
  width: auto;
  /* make sure input baseline aligns with prompt */
  vertical-align: baseline;
  /* padding + margin = 0.5em between prompt and cursor */
  padding: 0em 0.25em;
  margin: 0em 0.25em;
}
input.raw_input:focus {
  box-shadow: none;
}
p.p-space {
  margin-bottom: 10px;
}
div.output_unrecognized {
  padding: 5px;
  font-weight: bold;
  color: red;
}
div.output_unrecognized a {
  color: inherit;
  text-decoration: none;
}
div.output_unrecognized a:hover {
  color: inherit;
  text-decoration: none;
}
.rendered_html {
  color: #000;
  /* any extras will just be numbers: */
}
.rendered_html em {
  font-style: italic;
}
.rendered_html strong {
  font-weight: bold;
}
.rendered_html u {
  text-decoration: underline;
}
.rendered_html :link {
  text-decoration: underline;
}
.rendered_html :visited {
  text-decoration: underline;
}
.rendered_html h1 {
  font-size: 185.7%;
  margin: 1.08em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
}
.rendered_html h2 {
  font-size: 157.1%;
  margin: 1.27em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
}
.rendered_html h3 {
  font-size: 128.6%;
  margin: 1.55em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
}
.rendered_html h4 {
  font-size: 100%;
  margin: 2em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
}
.rendered_html h5 {
  font-size: 100%;
  margin: 2em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
  font-style: italic;
}
.rendered_html h6 {
  font-size: 100%;
  margin: 2em 0 0 0;
  font-weight: bold;
  line-height: 1.0;
  font-style: italic;
}
.rendered_html h1:first-child {
  margin-top: 0.538em;
}
.rendered_html h2:first-child {
  margin-top: 0.636em;
}
.rendered_html h3:first-child {
  margin-top: 0.777em;
}
.rendered_html h4:first-child {
  margin-top: 1em;
}
.rendered_html h5:first-child {
  margin-top: 1em;
}
.rendered_html h6:first-child {
  margin-top: 1em;
}
.rendered_html ul {
  list-style: disc;
  margin: 0em 2em;
  padding-left: 0px;
}
.rendered_html ul ul {
  list-style: square;
  margin: 0em 2em;
}
.rendered_html ul ul ul {
  list-style: circle;
  margin: 0em 2em;
}
.rendered_html ol {
  list-style: decimal;
  margin: 0em 2em;
  padding-left: 0px;
}
.rendered_html ol ol {
  list-style: upper-alpha;
  margin: 0em 2em;
}
.rendered_html ol ol ol {
  list-style: lower-alpha;
  margin: 0em 2em;
}
.rendered_html ol ol ol ol {
  list-style: lower-roman;
  margin: 0em 2em;
}
.rendered_html ol ol ol ol ol {
  list-style: decimal;
  margin: 0em 2em;
}
.rendered_html * + ul {
  margin-top: 1em;
}
.rendered_html * + ol {
  margin-top: 1em;
}
.rendered_html hr {
  color: black;
  background-color: black;
}
.rendered_html pre {
  margin: 1em 2em;
}
.rendered_html pre,
.rendered_html code {
  border: 0;
  background-color: #fff;
  color: #000;
  font-size: 100%;
  padding: 0px;
}
.rendered_html blockquote {
  margin: 1em 2em;
}
.rendered_html table {
  margin-left: auto;
  margin-right: auto;
  border: 1px solid black;
  border-collapse: collapse;
}
.rendered_html tr,
.rendered_html th,
.rendered_html td {
  border: 1px solid black;
  border-collapse: collapse;
  margin: 1em 2em;
}
.rendered_html td,
.rendered_html th {
  text-align: left;
  vertical-align: middle;
  padding: 4px;
}
.rendered_html th {
  font-weight: bold;
}
.rendered_html * + table {
  margin-top: 1em;
}
.rendered_html p {
  text-align: left;
}
.rendered_html * + p {
  margin-top: 1em;
}
.rendered_html img {
  display: block;
  margin-left: auto;
  margin-right: auto;
}
.rendered_html * + img {
  margin-top: 1em;
}
.rendered_html img,
.rendered_html svg {
  max-width: 100%;
  height: auto;
}
.rendered_html img.unconfined,
.rendered_html svg.unconfined {
  max-width: none;
}
div.text_cell {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
}
@media (max-width: 540px) {
  div.text_cell > div.prompt {
    display: none;
  }
}
div.text_cell_render {
  /*font-family: "Helvetica Neue", Arial, Helvetica, Geneva, sans-serif;*/
  outline: none;
  resize: none;
  width: inherit;
  border-style: none;
  padding: 0.5em 0.5em 0.5em 0.4em;
  color: #000;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}
a.anchor-link:link {
  text-decoration: none;
  padding: 0px 20px;
  visibility: hidden;
}
h1:hover .anchor-link,
h2:hover .anchor-link,
h3:hover .anchor-link,
h4:hover .anchor-link,
h5:hover .anchor-link,
h6:hover .anchor-link {
  visibility: visible;
}
.text_cell.rendered .input_area {
  display: none;
}
.text_cell.rendered .rendered_html {
  overflow-x: auto;
  overflow-y: hidden;
}
.text_cell.unrendered .text_cell_render {
  display: none;
}
.cm-header-1,
.cm-header-2,
.cm-header-3,
.cm-header-4,
.cm-header-5,
.cm-header-6 {
  font-weight: bold;
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
}
.cm-header-1 {
  font-size: 185.7%;
}
.cm-header-2 {
  font-size: 157.1%;
}
.cm-header-3 {
  font-size: 128.6%;
}
.cm-header-4 {
  font-size: 110%;
}
.cm-header-5 {
  font-size: 100%;
  font-style: italic;
}
.cm-header-6 {
  font-size: 100%;
  font-style: italic;
}
/*!
*
* IPython notebook webapp
*
*/
@media (max-width: 767px) {
  .notebook_app {
    padding-left: 0px;
    padding-right: 0px;
  }
}
#ipython-main-app {
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  height: 100%;
}
div#notebook_panel {
  margin: 0px;
  padding: 0px;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  height: 100%;
}
div#notebook {
  font-size: 14px;
  line-height: 20px;
  overflow-y: hidden;
  overflow-x: auto;
  width: 100%;
  /* This spaces the page away from the edge of the notebook area */
  padding-top: 20px;
  margin: 0px;
  outline: none;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  min-height: 100%;
}
@media not print {
  #notebook-container {
    padding: 15px;
    background-color: #fff;
    min-height: 0;
    -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
    box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  }
}
@media print {
  #notebook-container {
    width: 100%;
  }
}
div.ui-widget-content {
  border: 1px solid #ababab;
  outline: none;
}
pre.dialog {
  background-color: #f7f7f7;
  border: 1px solid #ddd;
  border-radius: 2px;
  padding: 0.4em;
  padding-left: 2em;
}
p.dialog {
  padding: 0.2em;
}
/* Word-wrap output correctly.  This is the CSS3 spelling, though Firefox seems
   to not honor it correctly.  Webkit browsers (Chrome, rekonq, Safari) do.
 */
pre,
code,
kbd,
samp {
  white-space: pre-wrap;
}
#fonttest {
  font-family: monospace;
}
p {
  margin-bottom: 0;
}
.end_space {
  min-height: 100px;
  transition: height .2s ease;
}
.notebook_app > #header {
  -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
}
@media not print {
  .notebook_app {
    background-color: #EEE;
  }
}
kbd {
  border-style: solid;
  border-width: 1px;
  box-shadow: none;
  margin: 2px;
  padding-left: 2px;
  padding-right: 2px;
  padding-top: 1px;
  padding-bottom: 1px;
}
/* CSS for the cell toolbar */
.celltoolbar {
  border: thin solid #CFCFCF;
  border-bottom: none;
  background: #EEE;
  border-radius: 2px 2px 0px 0px;
  width: 100%;
  height: 29px;
  padding-right: 4px;
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
  /* Old browsers */
  -webkit-box-pack: end;
  -moz-box-pack: end;
  box-pack: end;
  /* Modern browsers */
  justify-content: flex-end;
  display: -webkit-flex;
}
@media print {
  .celltoolbar {
    display: none;
  }
}
.ctb_hideshow {
  display: none;
  vertical-align: bottom;
}
/* ctb_show is added to the ctb_hideshow div to show the cell toolbar.
   Cell toolbars are only shown when the ctb_global_show class is also set.
*/
.ctb_global_show .ctb_show.ctb_hideshow {
  display: block;
}
.ctb_global_show .ctb_show + .input_area,
.ctb_global_show .ctb_show + div.text_cell_input,
.ctb_global_show .ctb_show ~ div.text_cell_render {
  border-top-right-radius: 0px;
  border-top-left-radius: 0px;
}
.ctb_global_show .ctb_show ~ div.text_cell_render {
  border: 1px solid #cfcfcf;
}
.celltoolbar {
  font-size: 87%;
  padding-top: 3px;
}
.celltoolbar select {
  display: block;
  width: 100%;
  height: 32px;
  padding: 6px 12px;
  font-size: 13px;
  line-height: 1.42857143;
  color: #555555;
  background-color: #fff;
  background-image: none;
  border: 1px solid #ccc;
  border-radius: 2px;
  -webkit-box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  box-shadow: inset 0 1px 1px rgba(0, 0, 0, 0.075);
  -webkit-transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
  -o-transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
  transition: border-color ease-in-out .15s, box-shadow ease-in-out .15s;
  height: 30px;
  padding: 5px 10px;
  font-size: 12px;
  line-height: 1.5;
  border-radius: 1px;
  width: inherit;
  font-size: inherit;
  height: 22px;
  padding: 0px;
  display: inline-block;
}
.celltoolbar select:focus {
  border-color: #66afe9;
  outline: 0;
  -webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.075), 0 0 8px rgba(102, 175, 233, 0.6);
  box-shadow: inset 0 1px 1px rgba(0,0,0,.075), 0 0 8px rgba(102, 175, 233, 0.6);
}
.celltoolbar select::-moz-placeholder {
  color: #999;
  opacity: 1;
}
.celltoolbar select:-ms-input-placeholder {
  color: #999;
}
.celltoolbar select::-webkit-input-placeholder {
  color: #999;
}
.celltoolbar select::-ms-expand {
  border: 0;
  background-color: transparent;
}
.celltoolbar select[disabled],
.celltoolbar select[readonly],
fieldset[disabled] .celltoolbar select {
  background-color: #eeeeee;
  opacity: 1;
}
.celltoolbar select[disabled],
fieldset[disabled] .celltoolbar select {
  cursor: not-allowed;
}
textarea.celltoolbar select {
  height: auto;
}
select.celltoolbar select {
  height: 30px;
  line-height: 30px;
}
textarea.celltoolbar select,
select[multiple].celltoolbar select {
  height: auto;
}
.celltoolbar label {
  margin-left: 5px;
  margin-right: 5px;
}
.completions {
  position: absolute;
  z-index: 110;
  overflow: hidden;
  border: 1px solid #ababab;
  border-radius: 2px;
  -webkit-box-shadow: 0px 6px 10px -1px #adadad;
  box-shadow: 0px 6px 10px -1px #adadad;
  line-height: 1;
}
.completions select {
  background: white;
  outline: none;
  border: none;
  padding: 0px;
  margin: 0px;
  overflow: auto;
  font-family: monospace;
  font-size: 110%;
  color: #000;
  width: auto;
}
.completions select option.context {
  color: #286090;
}
#kernel_logo_widget {
  float: right !important;
  float: right;
}
#kernel_logo_widget .current_kernel_logo {
  display: none;
  margin-top: -1px;
  margin-bottom: -1px;
  width: 32px;
  height: 32px;
}
#menubar {
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  margin-top: 1px;
}
#menubar .navbar {
  border-top: 1px;
  border-radius: 0px 0px 2px 2px;
  margin-bottom: 0px;
}
#menubar .navbar-toggle {
  float: left;
  padding-top: 7px;
  padding-bottom: 7px;
  border: none;
}
#menubar .navbar-collapse {
  clear: left;
}
.nav-wrapper {
  border-bottom: 1px solid #e7e7e7;
}
i.menu-icon {
  padding-top: 4px;
}
ul#help_menu li a {
  overflow: hidden;
  padding-right: 2.2em;
}
ul#help_menu li a i {
  margin-right: -1.2em;
}
.dropdown-submenu {
  position: relative;
}
.dropdown-submenu > .dropdown-menu {
  top: 0;
  left: 100%;
  margin-top: -6px;
  margin-left: -1px;
}
.dropdown-submenu:hover > .dropdown-menu {
  display: block;
}
.dropdown-submenu > a:after {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  display: block;
  content: "\f0da";
  float: right;
  color: #333333;
  margin-top: 2px;
  margin-right: -10px;
}
.dropdown-submenu > a:after.pull-left {
  margin-right: .3em;
}
.dropdown-submenu > a:after.pull-right {
  margin-left: .3em;
}
.dropdown-submenu:hover > a:after {
  color: #262626;
}
.dropdown-submenu.pull-left {
  float: none;
}
.dropdown-submenu.pull-left > .dropdown-menu {
  left: -100%;
  margin-left: 10px;
}
#notification_area {
  float: right !important;
  float: right;
  z-index: 10;
}
.indicator_area {
  float: right !important;
  float: right;
  color: #777;
  margin-left: 5px;
  margin-right: 5px;
  width: 11px;
  z-index: 10;
  text-align: center;
  width: auto;
}
#kernel_indicator {
  float: right !important;
  float: right;
  color: #777;
  margin-left: 5px;
  margin-right: 5px;
  width: 11px;
  z-index: 10;
  text-align: center;
  width: auto;
  border-left: 1px solid;
}
#kernel_indicator .kernel_indicator_name {
  padding-left: 5px;
  padding-right: 5px;
}
#modal_indicator {
  float: right !important;
  float: right;
  color: #777;
  margin-left: 5px;
  margin-right: 5px;
  width: 11px;
  z-index: 10;
  text-align: center;
  width: auto;
}
#readonly-indicator {
  float: right !important;
  float: right;
  color: #777;
  margin-left: 5px;
  margin-right: 5px;
  width: 11px;
  z-index: 10;
  text-align: center;
  width: auto;
  margin-top: 2px;
  margin-bottom: 0px;
  margin-left: 0px;
  margin-right: 0px;
  display: none;
}
.modal_indicator:before {
  width: 1.28571429em;
  text-align: center;
}
.edit_mode .modal_indicator:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f040";
}
.edit_mode .modal_indicator:before.pull-left {
  margin-right: .3em;
}
.edit_mode .modal_indicator:before.pull-right {
  margin-left: .3em;
}
.command_mode .modal_indicator:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: ' ';
}
.command_mode .modal_indicator:before.pull-left {
  margin-right: .3em;
}
.command_mode .modal_indicator:before.pull-right {
  margin-left: .3em;
}
.kernel_idle_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f10c";
}
.kernel_idle_icon:before.pull-left {
  margin-right: .3em;
}
.kernel_idle_icon:before.pull-right {
  margin-left: .3em;
}
.kernel_busy_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f111";
}
.kernel_busy_icon:before.pull-left {
  margin-right: .3em;
}
.kernel_busy_icon:before.pull-right {
  margin-left: .3em;
}
.kernel_dead_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f1e2";
}
.kernel_dead_icon:before.pull-left {
  margin-right: .3em;
}
.kernel_dead_icon:before.pull-right {
  margin-left: .3em;
}
.kernel_disconnected_icon:before {
  display: inline-block;
  font: normal normal normal 14px/1 FontAwesome;
  font-size: inherit;
  text-rendering: auto;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  content: "\f127";
}
.kernel_disconnected_icon:before.pull-left {
  margin-right: .3em;
}
.kernel_disconnected_icon:before.pull-right {
  margin-left: .3em;
}
.notification_widget {
  color: #777;
  z-index: 10;
  background: rgba(240, 240, 240, 0.5);
  margin-right: 4px;
  color: #333;
  background-color: #fff;
  border-color: #ccc;
}
.notification_widget:focus,
.notification_widget.focus {
  color: #333;
  background-color: #e6e6e6;
  border-color: #8c8c8c;
}
.notification_widget:hover {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
.notification_widget:active,
.notification_widget.active,
.open > .dropdown-toggle.notification_widget {
  color: #333;
  background-color: #e6e6e6;
  border-color: #adadad;
}
.notification_widget:active:hover,
.notification_widget.active:hover,
.open > .dropdown-toggle.notification_widget:hover,
.notification_widget:active:focus,
.notification_widget.active:focus,
.open > .dropdown-toggle.notification_widget:focus,
.notification_widget:active.focus,
.notification_widget.active.focus,
.open > .dropdown-toggle.notification_widget.focus {
  color: #333;
  background-color: #d4d4d4;
  border-color: #8c8c8c;
}
.notification_widget:active,
.notification_widget.active,
.open > .dropdown-toggle.notification_widget {
  background-image: none;
}
.notification_widget.disabled:hover,
.notification_widget[disabled]:hover,
fieldset[disabled] .notification_widget:hover,
.notification_widget.disabled:focus,
.notification_widget[disabled]:focus,
fieldset[disabled] .notification_widget:focus,
.notification_widget.disabled.focus,
.notification_widget[disabled].focus,
fieldset[disabled] .notification_widget.focus {
  background-color: #fff;
  border-color: #ccc;
}
.notification_widget .badge {
  color: #fff;
  background-color: #333;
}
.notification_widget.warning {
  color: #fff;
  background-color: #f0ad4e;
  border-color: #eea236;
}
.notification_widget.warning:focus,
.notification_widget.warning.focus {
  color: #fff;
  background-color: #ec971f;
  border-color: #985f0d;
}
.notification_widget.warning:hover {
  color: #fff;
  background-color: #ec971f;
  border-color: #d58512;
}
.notification_widget.warning:active,
.notification_widget.warning.active,
.open > .dropdown-toggle.notification_widget.warning {
  color: #fff;
  background-color: #ec971f;
  border-color: #d58512;
}
.notification_widget.warning:active:hover,
.notification_widget.warning.active:hover,
.open > .dropdown-toggle.notification_widget.warning:hover,
.notification_widget.warning:active:focus,
.notification_widget.warning.active:focus,
.open > .dropdown-toggle.notification_widget.warning:focus,
.notification_widget.warning:active.focus,
.notification_widget.warning.active.focus,
.open > .dropdown-toggle.notification_widget.warning.focus {
  color: #fff;
  background-color: #d58512;
  border-color: #985f0d;
}
.notification_widget.warning:active,
.notification_widget.warning.active,
.open > .dropdown-toggle.notification_widget.warning {
  background-image: none;
}
.notification_widget.warning.disabled:hover,
.notification_widget.warning[disabled]:hover,
fieldset[disabled] .notification_widget.warning:hover,
.notification_widget.warning.disabled:focus,
.notification_widget.warning[disabled]:focus,
fieldset[disabled] .notification_widget.warning:focus,
.notification_widget.warning.disabled.focus,
.notification_widget.warning[disabled].focus,
fieldset[disabled] .notification_widget.warning.focus {
  background-color: #f0ad4e;
  border-color: #eea236;
}
.notification_widget.warning .badge {
  color: #f0ad4e;
  background-color: #fff;
}
.notification_widget.success {
  color: #fff;
  background-color: #5cb85c;
  border-color: #4cae4c;
}
.notification_widget.success:focus,
.notification_widget.success.focus {
  color: #fff;
  background-color: #449d44;
  border-color: #255625;
}
.notification_widget.success:hover {
  color: #fff;
  background-color: #449d44;
  border-color: #398439;
}
.notification_widget.success:active,
.notification_widget.success.active,
.open > .dropdown-toggle.notification_widget.success {
  color: #fff;
  background-color: #449d44;
  border-color: #398439;
}
.notification_widget.success:active:hover,
.notification_widget.success.active:hover,
.open > .dropdown-toggle.notification_widget.success:hover,
.notification_widget.success:active:focus,
.notification_widget.success.active:focus,
.open > .dropdown-toggle.notification_widget.success:focus,
.notification_widget.success:active.focus,
.notification_widget.success.active.focus,
.open > .dropdown-toggle.notification_widget.success.focus {
  color: #fff;
  background-color: #398439;
  border-color: #255625;
}
.notification_widget.success:active,
.notification_widget.success.active,
.open > .dropdown-toggle.notification_widget.success {
  background-image: none;
}
.notification_widget.success.disabled:hover,
.notification_widget.success[disabled]:hover,
fieldset[disabled] .notification_widget.success:hover,
.notification_widget.success.disabled:focus,
.notification_widget.success[disabled]:focus,
fieldset[disabled] .notification_widget.success:focus,
.notification_widget.success.disabled.focus,
.notification_widget.success[disabled].focus,
fieldset[disabled] .notification_widget.success.focus {
  background-color: #5cb85c;
  border-color: #4cae4c;
}
.notification_widget.success .badge {
  color: #5cb85c;
  background-color: #fff;
}
.notification_widget.info {
  color: #fff;
  background-color: #5bc0de;
  border-color: #46b8da;
}
.notification_widget.info:focus,
.notification_widget.info.focus {
  color: #fff;
  background-color: #31b0d5;
  border-color: #1b6d85;
}
.notification_widget.info:hover {
  color: #fff;
  background-color: #31b0d5;
  border-color: #269abc;
}
.notification_widget.info:active,
.notification_widget.info.active,
.open > .dropdown-toggle.notification_widget.info {
  color: #fff;
  background-color: #31b0d5;
  border-color: #269abc;
}
.notification_widget.info:active:hover,
.notification_widget.info.active:hover,
.open > .dropdown-toggle.notification_widget.info:hover,
.notification_widget.info:active:focus,
.notification_widget.info.active:focus,
.open > .dropdown-toggle.notification_widget.info:focus,
.notification_widget.info:active.focus,
.notification_widget.info.active.focus,
.open > .dropdown-toggle.notification_widget.info.focus {
  color: #fff;
  background-color: #269abc;
  border-color: #1b6d85;
}
.notification_widget.info:active,
.notification_widget.info.active,
.open > .dropdown-toggle.notification_widget.info {
  background-image: none;
}
.notification_widget.info.disabled:hover,
.notification_widget.info[disabled]:hover,
fieldset[disabled] .notification_widget.info:hover,
.notification_widget.info.disabled:focus,
.notification_widget.info[disabled]:focus,
fieldset[disabled] .notification_widget.info:focus,
.notification_widget.info.disabled.focus,
.notification_widget.info[disabled].focus,
fieldset[disabled] .notification_widget.info.focus {
  background-color: #5bc0de;
  border-color: #46b8da;
}
.notification_widget.info .badge {
  color: #5bc0de;
  background-color: #fff;
}
.notification_widget.danger {
  color: #fff;
  background-color: #d9534f;
  border-color: #d43f3a;
}
.notification_widget.danger:focus,
.notification_widget.danger.focus {
  color: #fff;
  background-color: #c9302c;
  border-color: #761c19;
}
.notification_widget.danger:hover {
  color: #fff;
  background-color: #c9302c;
  border-color: #ac2925;
}
.notification_widget.danger:active,
.notification_widget.danger.active,
.open > .dropdown-toggle.notification_widget.danger {
  color: #fff;
  background-color: #c9302c;
  border-color: #ac2925;
}
.notification_widget.danger:active:hover,
.notification_widget.danger.active:hover,
.open > .dropdown-toggle.notification_widget.danger:hover,
.notification_widget.danger:active:focus,
.notification_widget.danger.active:focus,
.open > .dropdown-toggle.notification_widget.danger:focus,
.notification_widget.danger:active.focus,
.notification_widget.danger.active.focus,
.open > .dropdown-toggle.notification_widget.danger.focus {
  color: #fff;
  background-color: #ac2925;
  border-color: #761c19;
}
.notification_widget.danger:active,
.notification_widget.danger.active,
.open > .dropdown-toggle.notification_widget.danger {
  background-image: none;
}
.notification_widget.danger.disabled:hover,
.notification_widget.danger[disabled]:hover,
fieldset[disabled] .notification_widget.danger:hover,
.notification_widget.danger.disabled:focus,
.notification_widget.danger[disabled]:focus,
fieldset[disabled] .notification_widget.danger:focus,
.notification_widget.danger.disabled.focus,
.notification_widget.danger[disabled].focus,
fieldset[disabled] .notification_widget.danger.focus {
  background-color: #d9534f;
  border-color: #d43f3a;
}
.notification_widget.danger .badge {
  color: #d9534f;
  background-color: #fff;
}
div#pager {
  background-color: #fff;
  font-size: 14px;
  line-height: 20px;
  overflow: hidden;
  display: none;
  position: fixed;
  bottom: 0px;
  width: 100%;
  max-height: 50%;
  padding-top: 8px;
  -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  /* Display over codemirror */
  z-index: 100;
  /* Hack which prevents jquery ui resizable from changing top. */
  top: auto !important;
}
div#pager pre {
  line-height: 1.21429em;
  color: #000;
  background-color: #f7f7f7;
  padding: 0.4em;
}
div#pager #pager-button-area {
  position: absolute;
  top: 8px;
  right: 20px;
}
div#pager #pager-contents {
  position: relative;
  overflow: auto;
  width: 100%;
  height: 100%;
}
div#pager #pager-contents #pager-container {
  position: relative;
  padding: 15px 0px;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}
div#pager .ui-resizable-handle {
  top: 0px;
  height: 8px;
  background: #f7f7f7;
  border-top: 1px solid #cfcfcf;
  border-bottom: 1px solid #cfcfcf;
  /* This injects handle bars (a short, wide = symbol) for 
        the resize handle. */
}
div#pager .ui-resizable-handle::after {
  content: '';
  top: 2px;
  left: 50%;
  height: 3px;
  width: 30px;
  margin-left: -15px;
  position: absolute;
  border-top: 1px solid #cfcfcf;
}
.quickhelp {
  /* Old browsers */
  display: -webkit-box;
  -webkit-box-orient: horizontal;
  -webkit-box-align: stretch;
  display: -moz-box;
  -moz-box-orient: horizontal;
  -moz-box-align: stretch;
  display: box;
  box-orient: horizontal;
  box-align: stretch;
  /* Modern browsers */
  display: flex;
  flex-direction: row;
  align-items: stretch;
  line-height: 1.8em;
}
.shortcut_key {
  display: inline-block;
  width: 21ex;
  text-align: right;
  font-family: monospace;
}
.shortcut_descr {
  display: inline-block;
  /* Old browsers */
  -webkit-box-flex: 1;
  -moz-box-flex: 1;
  box-flex: 1;
  /* Modern browsers */
  flex: 1;
}
span.save_widget {
  margin-top: 6px;
}
span.save_widget span.filename {
  height: 1em;
  line-height: 1em;
  padding: 3px;
  margin-left: 16px;
  border: none;
  font-size: 146.5%;
  border-radius: 2px;
}
span.save_widget span.filename:hover {
  background-color: #e6e6e6;
}
span.checkpoint_status,
span.autosave_status {
  font-size: small;
}
@media (max-width: 767px) {
  span.save_widget {
    font-size: small;
  }
  span.checkpoint_status,
  span.autosave_status {
    display: none;
  }
}
@media (min-width: 768px) and (max-width: 991px) {
  span.checkpoint_status {
    display: none;
  }
  span.autosave_status {
    font-size: x-small;
  }
}
.toolbar {
  padding: 0px;
  margin-left: -5px;
  margin-top: 2px;
  margin-bottom: 5px;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}
.toolbar select,
.toolbar label {
  width: auto;
  vertical-align: middle;
  margin-right: 2px;
  margin-bottom: 0px;
  display: inline;
  font-size: 92%;
  margin-left: 0.3em;
  margin-right: 0.3em;
  padding: 0px;
  padding-top: 3px;
}
.toolbar .btn {
  padding: 2px 8px;
}
.toolbar .btn-group {
  margin-top: 0px;
  margin-left: 5px;
}
#maintoolbar {
  margin-bottom: -3px;
  margin-top: -8px;
  border: 0px;
  min-height: 27px;
  margin-left: 0px;
  padding-top: 11px;
  padding-bottom: 3px;
}
#maintoolbar .navbar-text {
  float: none;
  vertical-align: middle;
  text-align: right;
  margin-left: 5px;
  margin-right: 0px;
  margin-top: 0px;
}
.select-xs {
  height: 24px;
}
.pulse,
.dropdown-menu > li > a.pulse,
li.pulse > a.dropdown-toggle,
li.pulse.open > a.dropdown-toggle {
  background-color: #F37626;
  color: white;
}
/**
 * Primary styles
 *
 * Author: Jupyter Development Team
 */
/** WARNING IF YOU ARE EDITTING THIS FILE, if this is a .css file, It has a lot
 * of chance of beeing generated from the ../less/[samename].less file, you can
 * try to get back the less file by reverting somme commit in history
 **/
/*
 * We'll try to get something pretty, so we
 * have some strange css to have the scroll bar on
 * the left with fix button on the top right of the tooltip
 */
@-moz-keyframes fadeOut {
  from {
    opacity: 1;
  }
  to {
    opacity: 0;
  }
}
@-webkit-keyframes fadeOut {
  from {
    opacity: 1;
  }
  to {
    opacity: 0;
  }
}
@-moz-keyframes fadeIn {
  from {
    opacity: 0;
  }
  to {
    opacity: 1;
  }
}
@-webkit-keyframes fadeIn {
  from {
    opacity: 0;
  }
  to {
    opacity: 1;
  }
}
/*properties of tooltip after "expand"*/
.bigtooltip {
  overflow: auto;
  height: 200px;
  -webkit-transition-property: height;
  -webkit-transition-duration: 500ms;
  -moz-transition-property: height;
  -moz-transition-duration: 500ms;
  transition-property: height;
  transition-duration: 500ms;
}
/*properties of tooltip before "expand"*/
.smalltooltip {
  -webkit-transition-property: height;
  -webkit-transition-duration: 500ms;
  -moz-transition-property: height;
  -moz-transition-duration: 500ms;
  transition-property: height;
  transition-duration: 500ms;
  text-overflow: ellipsis;
  overflow: hidden;
  height: 80px;
}
.tooltipbuttons {
  position: absolute;
  padding-right: 15px;
  top: 0px;
  right: 0px;
}
.tooltiptext {
  /*avoid the button to overlap on some docstring*/
  padding-right: 30px;
}
.ipython_tooltip {
  max-width: 700px;
  /*fade-in animation when inserted*/
  -webkit-animation: fadeOut 400ms;
  -moz-animation: fadeOut 400ms;
  animation: fadeOut 400ms;
  -webkit-animation: fadeIn 400ms;
  -moz-animation: fadeIn 400ms;
  animation: fadeIn 400ms;
  vertical-align: middle;
  background-color: #f7f7f7;
  overflow: visible;
  border: #ababab 1px solid;
  outline: none;
  padding: 3px;
  margin: 0px;
  padding-left: 7px;
  font-family: monospace;
  min-height: 50px;
  -moz-box-shadow: 0px 6px 10px -1px #adadad;
  -webkit-box-shadow: 0px 6px 10px -1px #adadad;
  box-shadow: 0px 6px 10px -1px #adadad;
  border-radius: 2px;
  position: absolute;
  z-index: 1000;
}
.ipython_tooltip a {
  float: right;
}
.ipython_tooltip .tooltiptext pre {
  border: 0;
  border-radius: 0;
  font-size: 100%;
  background-color: #f7f7f7;
}
.pretooltiparrow {
  left: 0px;
  margin: 0px;
  top: -16px;
  width: 40px;
  height: 16px;
  overflow: hidden;
  position: absolute;
}
.pretooltiparrow:before {
  background-color: #f7f7f7;
  border: 1px #ababab solid;
  z-index: 11;
  content: "";
  position: absolute;
  left: 15px;
  top: 10px;
  width: 25px;
  height: 25px;
  -webkit-transform: rotate(45deg);
  -moz-transform: rotate(45deg);
  -ms-transform: rotate(45deg);
  -o-transform: rotate(45deg);
}
ul.typeahead-list i {
  margin-left: -10px;
  width: 18px;
}
ul.typeahead-list {
  max-height: 80vh;
  overflow: auto;
}
ul.typeahead-list > li > a {
  /** Firefox bug **/
  /* see https://github.com/jupyter/notebook/issues/559 */
  white-space: normal;
}
.cmd-palette .modal-body {
  padding: 7px;
}
.cmd-palette form {
  background: white;
}
.cmd-palette input {
  outline: none;
}
.no-shortcut {
  display: none;
}
.command-shortcut:before {
  content: "(command)";
  padding-right: 3px;
  color: #777777;
}
.edit-shortcut:before {
  content: "(edit)";
  padding-right: 3px;
  color: #777777;
}
#find-and-replace #replace-preview .match,
#find-and-replace #replace-preview .insert {
  background-color: #BBDEFB;
  border-color: #90CAF9;
  border-style: solid;
  border-width: 1px;
  border-radius: 0px;
}
#find-and-replace #replace-preview .replace .match {
  background-color: #FFCDD2;
  border-color: #EF9A9A;
  border-radius: 0px;
}
#find-and-replace #replace-preview .replace .insert {
  background-color: #C8E6C9;
  border-color: #A5D6A7;
  border-radius: 0px;
}
#find-and-replace #replace-preview {
  max-height: 60vh;
  overflow: auto;
}
#find-and-replace #replace-preview pre {
  padding: 5px 10px;
}
.terminal-app {
  background: #EEE;
}
.terminal-app #header {
  background: #fff;
  -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
  box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.2);
}
.terminal-app .terminal {
  width: 100%;
  float: left;
  font-family: monospace;
  color: white;
  background: black;
  padding: 0.4em;
  border-radius: 2px;
  -webkit-box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.4);
  box-shadow: 0px 0px 12px 1px rgba(87, 87, 87, 0.4);
}
.terminal-app .terminal,
.terminal-app .terminal dummy-screen {
  line-height: 1em;
  font-size: 14px;
}
.terminal-app .terminal .xterm-rows {
  padding: 10px;
}
.terminal-app .terminal-cursor {
  color: black;
  background: white;
}
.terminal-app #terminado-container {
  margin-top: 20px;
}
/*# sourceMappingURL=style.min.css.map */
    </style>
<style type="text/css">
    .highlight .hll { background-color: #ffffcc }
.highlight  { background: #f8f8f8; }
.highlight .c { color: #408080; font-style: italic } /* Comment */
.highlight .err { border: 1px solid #FF0000 } /* Error */
.highlight .k { color: #008000; font-weight: bold } /* Keyword */
.highlight .o { color: #666666 } /* Operator */
.highlight .ch { color: #408080; font-style: italic } /* Comment.Hashbang */
.highlight .cm { color: #408080; font-style: italic } /* Comment.Multiline */
.highlight .cp { color: #BC7A00 } /* Comment.Preproc */
.highlight .cpf { color: #408080; font-style: italic } /* Comment.PreprocFile */
.highlight .c1 { color: #408080; font-style: italic } /* Comment.Single */
.highlight .cs { color: #408080; font-style: italic } /* Comment.Special */
.highlight .gd { color: #A00000 } /* Generic.Deleted */
.highlight .ge { font-style: italic } /* Generic.Emph */
.highlight .gr { color: #FF0000 } /* Generic.Error */
.highlight .gh { color: #000080; font-weight: bold } /* Generic.Heading */
.highlight .gi { color: #00A000 } /* Generic.Inserted */
.highlight .go { color: #888888 } /* Generic.Output */
.highlight .gp { color: #000080; font-weight: bold } /* Generic.Prompt */
.highlight .gs { font-weight: bold } /* Generic.Strong */
.highlight .gu { color: #800080; font-weight: bold } /* Generic.Subheading */
.highlight .gt { color: #0044DD } /* Generic.Traceback */
.highlight .kc { color: #008000; font-weight: bold } /* Keyword.Constant */
.highlight .kd { color: #008000; font-weight: bold } /* Keyword.Declaration */
.highlight .kn { color: #008000; font-weight: bold } /* Keyword.Namespace */
.highlight .kp { color: #008000 } /* Keyword.Pseudo */
.highlight .kr { color: #008000; font-weight: bold } /* Keyword.Reserved */
.highlight .kt { color: #B00040 } /* Keyword.Type */
.highlight .m { color: #666666 } /* Literal.Number */
.highlight .s { color: #BA2121 } /* Literal.String */
.highlight .na { color: #7D9029 } /* Name.Attribute */
.highlight .nb { color: #008000 } /* Name.Builtin */
.highlight .nc { color: #0000FF; font-weight: bold } /* Name.Class */
.highlight .no { color: #880000 } /* Name.Constant */
.highlight .nd { color: #AA22FF } /* Name.Decorator */
.highlight .ni { color: #999999; font-weight: bold } /* Name.Entity */
.highlight .ne { color: #D2413A; font-weight: bold } /* Name.Exception */
.highlight .nf { color: #0000FF } /* Name.Function */
.highlight .nl { color: #A0A000 } /* Name.Label */
.highlight .nn { color: #0000FF; font-weight: bold } /* Name.Namespace */
.highlight .nt { color: #008000; font-weight: bold } /* Name.Tag */
.highlight .nv { color: #19177C } /* Name.Variable */
.highlight .ow { color: #AA22FF; font-weight: bold } /* Operator.Word */
.highlight .w { color: #bbbbbb } /* Text.Whitespace */
.highlight .mb { color: #666666 } /* Literal.Number.Bin */
.highlight .mf { color: #666666 } /* Literal.Number.Float */
.highlight .mh { color: #666666 } /* Literal.Number.Hex */
.highlight .mi { color: #666666 } /* Literal.Number.Integer */
.highlight .mo { color: #666666 } /* Literal.Number.Oct */
.highlight .sa { color: #BA2121 } /* Literal.String.Affix */
.highlight .sb { color: #BA2121 } /* Literal.String.Backtick */
.highlight .sc { color: #BA2121 } /* Literal.String.Char */
.highlight .dl { color: #BA2121 } /* Literal.String.Delimiter */
.highlight .sd { color: #BA2121; font-style: italic } /* Literal.String.Doc */
.highlight .s2 { color: #BA2121 } /* Literal.String.Double */
.highlight .se { color: #BB6622; font-weight: bold } /* Literal.String.Escape */
.highlight .sh { color: #BA2121 } /* Literal.String.Heredoc */
.highlight .si { color: #BB6688; font-weight: bold } /* Literal.String.Interpol */
.highlight .sx { color: #008000 } /* Literal.String.Other */
.highlight .sr { color: #BB6688 } /* Literal.String.Regex */
.highlight .s1 { color: #BA2121 } /* Literal.String.Single */
.highlight .ss { color: #19177C } /* Literal.String.Symbol */
.highlight .bp { color: #008000 } /* Name.Builtin.Pseudo */
.highlight .fm { color: #0000FF } /* Name.Function.Magic */
.highlight .vc { color: #19177C } /* Name.Variable.Class */
.highlight .vg { color: #19177C } /* Name.Variable.Global */
.highlight .vi { color: #19177C } /* Name.Variable.Instance */
.highlight .vm { color: #19177C } /* Name.Variable.Magic */
.highlight .il { color: #666666 } /* Literal.Number.Integer.Long */
    </style>
<style type="text/css">
    
/* Temporary definitions which will become obsolete with Notebook release 5.0 */
.ansi-black-fg { color: #3E424D; }
.ansi-black-bg { background-color: #3E424D; }
.ansi-black-intense-fg { color: #282C36; }
.ansi-black-intense-bg { background-color: #282C36; }
.ansi-red-fg { color: #E75C58; }
.ansi-red-bg { background-color: #E75C58; }
.ansi-red-intense-fg { color: #B22B31; }
.ansi-red-intense-bg { background-color: #B22B31; }
.ansi-green-fg { color: #00A250; }
.ansi-green-bg { background-color: #00A250; }
.ansi-green-intense-fg { color: #007427; }
.ansi-green-intense-bg { background-color: #007427; }
.ansi-yellow-fg { color: #DDB62B; }
.ansi-yellow-bg { background-color: #DDB62B; }
.ansi-yellow-intense-fg { color: #B27D12; }
.ansi-yellow-intense-bg { background-color: #B27D12; }
.ansi-blue-fg { color: #208FFB; }
.ansi-blue-bg { background-color: #208FFB; }
.ansi-blue-intense-fg { color: #0065CA; }
.ansi-blue-intense-bg { background-color: #0065CA; }
.ansi-magenta-fg { color: #D160C4; }
.ansi-magenta-bg { background-color: #D160C4; }
.ansi-magenta-intense-fg { color: #A03196; }
.ansi-magenta-intense-bg { background-color: #A03196; }
.ansi-cyan-fg { color: #60C6C8; }
.ansi-cyan-bg { background-color: #60C6C8; }
.ansi-cyan-intense-fg { color: #258F8F; }
.ansi-cyan-intense-bg { background-color: #258F8F; }
.ansi-white-fg { color: #C5C1B4; }
.ansi-white-bg { background-color: #C5C1B4; }
.ansi-white-intense-fg { color: #A1A6B2; }
.ansi-white-intense-bg { background-color: #A1A6B2; }

.ansi-bold { font-weight: bold; }

    </style>


<style type="text/css">
/* Overrides of notebook CSS for static HTML export */
body {
  overflow: visible;
  padding: 8px;
}

div#notebook {
  overflow: visible;
  border-top: none;
}@media print {
  div.cell {
    display: block;
    page-break-inside: avoid;
  } 
  div.output_wrapper { 
    display: block;
    page-break-inside: avoid; 
  }
  div.output { 
    display: block;
    page-break-inside: avoid; 
  }
}
</style>

<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="./mdm_model_formula_check_files/custom.css">

<!-- Loading mathjax macro -->
<!-- Load mathjax -->
    <script src="./mdm_model_formula_check_files/MathJax.js.download"></script>
    <!-- MathJax configuration -->
    <script type="text/x-mathjax-config;executed=true">
    MathJax.Hub.Config({
        tex2jax: {
            inlineMath: [ ['$','$'], ["\\(","\\)"] ],
            displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
            processEscapes: true,
            processEnvironments: true
        },
        // Center justify equations in code and markdown cells. Elsewhere
        // we use CSS to left justify single line equations in code cells.
        displayAlign: 'center',
        "HTML-CSS": {
            styles: {'.MathJax_Display': {"margin": 0}},
            linebreaks: { automatic: true }
        }
    });
    </script>
    <!-- End of mathjax configuration --><style type="text/css">.MathJax_Hover_Frame {border-radius: .25em; -webkit-border-radius: .25em; -moz-border-radius: .25em; -khtml-border-radius: .25em; box-shadow: 0px 0px 15px #83A; -webkit-box-shadow: 0px 0px 15px #83A; -moz-box-shadow: 0px 0px 15px #83A; -khtml-box-shadow: 0px 0px 15px #83A; border: 1px solid #A6D ! important; display: inline-block; position: absolute}
.MathJax_Menu_Button .MathJax_Hover_Arrow {position: absolute; cursor: pointer; display: inline-block; border: 2px solid #AAA; border-radius: 4px; -webkit-border-radius: 4px; -moz-border-radius: 4px; -khtml-border-radius: 4px; font-family: 'Courier New',Courier; font-size: 9px; color: #F0F0F0}
.MathJax_Menu_Button .MathJax_Hover_Arrow span {display: block; background-color: #AAA; border: 1px solid; border-radius: 3px; line-height: 0; padding: 4px}
.MathJax_Hover_Arrow:hover {color: white!important; border: 2px solid #CCC!important}
.MathJax_Hover_Arrow:hover span {background-color: #CCC!important}
</style><style type="text/css">#MathJax_About {position: fixed; left: 50%; width: auto; text-align: center; border: 3px outset; padding: 1em 2em; background-color: #DDDDDD; color: black; cursor: default; font-family: message-box; font-size: 120%; font-style: normal; text-indent: 0; text-transform: none; line-height: normal; letter-spacing: normal; word-spacing: normal; word-wrap: normal; white-space: nowrap; float: none; z-index: 201; border-radius: 15px; -webkit-border-radius: 15px; -moz-border-radius: 15px; -khtml-border-radius: 15px; box-shadow: 0px 10px 20px #808080; -webkit-box-shadow: 0px 10px 20px #808080; -moz-box-shadow: 0px 10px 20px #808080; -khtml-box-shadow: 0px 10px 20px #808080; filter: progid:DXImageTransform.Microsoft.dropshadow(OffX=2, OffY=2, Color='gray', Positive='true')}
#MathJax_About.MathJax_MousePost {outline: none}
.MathJax_Menu {position: absolute; background-color: white; color: black; width: auto; padding: 2px; border: 1px solid #CCCCCC; margin: 0; cursor: default; font: menu; text-align: left; text-indent: 0; text-transform: none; line-height: normal; letter-spacing: normal; word-spacing: normal; word-wrap: normal; white-space: nowrap; float: none; z-index: 201; box-shadow: 0px 10px 20px #808080; -webkit-box-shadow: 0px 10px 20px #808080; -moz-box-shadow: 0px 10px 20px #808080; -khtml-box-shadow: 0px 10px 20px #808080; filter: progid:DXImageTransform.Microsoft.dropshadow(OffX=2, OffY=2, Color='gray', Positive='true')}
.MathJax_MenuItem {padding: 2px 2em; background: transparent}
.MathJax_MenuArrow {position: absolute; right: .5em; padding-top: .25em; color: #666666; font-size: .75em}
.MathJax_MenuActive .MathJax_MenuArrow {color: white}
.MathJax_MenuArrow.RTL {left: .5em; right: auto}
.MathJax_MenuCheck {position: absolute; left: .7em}
.MathJax_MenuCheck.RTL {right: .7em; left: auto}
.MathJax_MenuRadioCheck {position: absolute; left: 1em}
.MathJax_MenuRadioCheck.RTL {right: 1em; left: auto}
.MathJax_MenuLabel {padding: 2px 2em 4px 1.33em; font-style: italic}
.MathJax_MenuRule {border-top: 1px solid #CCCCCC; margin: 4px 1px 0px}
.MathJax_MenuDisabled {color: GrayText}
.MathJax_MenuActive {background-color: Highlight; color: HighlightText}
.MathJax_MenuDisabled:focus, .MathJax_MenuLabel:focus {background-color: #E8E8E8}
.MathJax_ContextMenu:focus {outline: none}
.MathJax_ContextMenu .MathJax_MenuItem:focus {outline: none}
#MathJax_AboutClose {top: .2em; right: .2em}
.MathJax_Menu .MathJax_MenuClose {top: -10px; left: -10px}
.MathJax_MenuClose {position: absolute; cursor: pointer; display: inline-block; border: 2px solid #AAA; border-radius: 18px; -webkit-border-radius: 18px; -moz-border-radius: 18px; -khtml-border-radius: 18px; font-family: 'Courier New',Courier; font-size: 24px; color: #F0F0F0}
.MathJax_MenuClose span {display: block; background-color: #AAA; border: 1.5px solid; border-radius: 18px; -webkit-border-radius: 18px; -moz-border-radius: 18px; -khtml-border-radius: 18px; line-height: 0; padding: 8px 0 6px}
.MathJax_MenuClose:hover {color: white!important; border: 2px solid #CCC!important}
.MathJax_MenuClose:hover span {background-color: #CCC!important}
.MathJax_MenuClose:hover:focus {outline: none}
</style><style type="text/css">.MathJax_Preview .MJXf-math {color: inherit!important}
</style><style type="text/css">.MJX_Assistive_MathML {position: absolute!important; top: 0; left: 0; clip: rect(1px, 1px, 1px, 1px); padding: 1px 0 0 0!important; border: 0!important; height: 1px!important; width: 1px!important; overflow: hidden!important; display: block!important; -webkit-touch-callout: none; -webkit-user-select: none; -khtml-user-select: none; -moz-user-select: none; -ms-user-select: none; user-select: none}
.MJX_Assistive_MathML.MJX_Assistive_MathML_Block {width: 100%!important}
</style><style type="text/css">#MathJax_Zoom {position: absolute; background-color: #F0F0F0; overflow: auto; display: block; z-index: 301; padding: .5em; border: 1px solid black; margin: 0; font-weight: normal; font-style: normal; text-align: left; text-indent: 0; text-transform: none; line-height: normal; letter-spacing: normal; word-spacing: normal; word-wrap: normal; white-space: nowrap; float: none; -webkit-box-sizing: content-box; -moz-box-sizing: content-box; box-sizing: content-box; box-shadow: 5px 5px 15px #AAAAAA; -webkit-box-shadow: 5px 5px 15px #AAAAAA; -moz-box-shadow: 5px 5px 15px #AAAAAA; -khtml-box-shadow: 5px 5px 15px #AAAAAA; filter: progid:DXImageTransform.Microsoft.dropshadow(OffX=2, OffY=2, Color='gray', Positive='true')}
#MathJax_ZoomOverlay {position: absolute; left: 0; top: 0; z-index: 300; display: inline-block; width: 100%; height: 100%; border: 0; padding: 0; margin: 0; background-color: white; opacity: 0; filter: alpha(opacity=0)}
#MathJax_ZoomFrame {position: relative; display: inline-block; height: 0; width: 0}
#MathJax_ZoomEventTrap {position: absolute; left: 0; top: 0; z-index: 302; display: inline-block; border: 0; padding: 0; margin: 0; background-color: white; opacity: 0; filter: alpha(opacity=0)}
</style><style type="text/css">.MathJax_Preview {color: #888}
#MathJax_Message {position: fixed; left: 1em; bottom: 1.5em; background-color: #E6E6E6; border: 1px solid #959595; margin: 0px; padding: 2px 8px; z-index: 102; color: black; font-size: 80%; width: auto; white-space: nowrap}
#MathJax_MSIE_Frame {position: absolute; top: 0; left: 0; width: 0px; z-index: 101; border: 0px; margin: 0px; padding: 0px}
.MathJax_Error {color: #CC0000; font-style: italic}
</style><style type="text/css">.MJXp-script {font-size: .8em}
.MJXp-right {-webkit-transform-origin: right; -moz-transform-origin: right; -ms-transform-origin: right; -o-transform-origin: right; transform-origin: right}
.MJXp-bold {font-weight: bold}
.MJXp-italic {font-style: italic}
.MJXp-scr {font-family: MathJax_Script,'Times New Roman',Times,STIXGeneral,serif}
.MJXp-frak {font-family: MathJax_Fraktur,'Times New Roman',Times,STIXGeneral,serif}
.MJXp-sf {font-family: MathJax_SansSerif,'Times New Roman',Times,STIXGeneral,serif}
.MJXp-cal {font-family: MathJax_Caligraphic,'Times New Roman',Times,STIXGeneral,serif}
.MJXp-mono {font-family: MathJax_Typewriter,'Times New Roman',Times,STIXGeneral,serif}
.MJXp-largeop {font-size: 150%}
.MJXp-largeop.MJXp-int {vertical-align: -.2em}
.MJXp-math {display: inline-block; line-height: 1.2; text-indent: 0; font-family: 'Times New Roman',Times,STIXGeneral,serif; white-space: nowrap; border-collapse: collapse}
.MJXp-display {display: block; text-align: center; margin: 1em 0}
.MJXp-math span {display: inline-block}
.MJXp-box {display: block!important; text-align: center}
.MJXp-box:after {content: " "}
.MJXp-rule {display: block!important; margin-top: .1em}
.MJXp-char {display: block!important}
.MJXp-mo {margin: 0 .15em}
.MJXp-mfrac {margin: 0 .125em; vertical-align: .25em}
.MJXp-denom {display: inline-table!important; width: 100%}
.MJXp-denom > * {display: table-row!important}
.MJXp-surd {vertical-align: top}
.MJXp-surd > * {display: block!important}
.MJXp-script-box > *  {display: table!important; height: 50%}
.MJXp-script-box > * > * {display: table-cell!important; vertical-align: top}
.MJXp-script-box > *:last-child > * {vertical-align: bottom}
.MJXp-script-box > * > * > * {display: block!important}
.MJXp-mphantom {visibility: hidden}
.MJXp-munderover {display: inline-table!important}
.MJXp-over {display: inline-block!important; text-align: center}
.MJXp-over > * {display: block!important}
.MJXp-munderover > * {display: table-row!important}
.MJXp-mtable {vertical-align: .25em; margin: 0 .125em}
.MJXp-mtable > * {display: inline-table!important; vertical-align: middle}
.MJXp-mtr {display: table-row!important}
.MJXp-mtd {display: table-cell!important; text-align: center; padding: .5em 0 0 .5em}
.MJXp-mtr > .MJXp-mtd:first-child {padding-left: 0}
.MJXp-mtr:first-child > .MJXp-mtd {padding-top: 0}
.MJXp-mlabeledtr {display: table-row!important}
.MJXp-mlabeledtr > .MJXp-mtd:first-child {padding-left: 0}
.MJXp-mlabeledtr:first-child > .MJXp-mtd {padding-top: 0}
.MJXp-merror {background-color: #FFFF88; color: #CC0000; border: 1px solid #CC0000; padding: 1px 3px; font-style: normal; font-size: 90%}
.MJXp-scale0 {-webkit-transform: scaleX(.0); -moz-transform: scaleX(.0); -ms-transform: scaleX(.0); -o-transform: scaleX(.0); transform: scaleX(.0)}
.MJXp-scale1 {-webkit-transform: scaleX(.1); -moz-transform: scaleX(.1); -ms-transform: scaleX(.1); -o-transform: scaleX(.1); transform: scaleX(.1)}
.MJXp-scale2 {-webkit-transform: scaleX(.2); -moz-transform: scaleX(.2); -ms-transform: scaleX(.2); -o-transform: scaleX(.2); transform: scaleX(.2)}
.MJXp-scale3 {-webkit-transform: scaleX(.3); -moz-transform: scaleX(.3); -ms-transform: scaleX(.3); -o-transform: scaleX(.3); transform: scaleX(.3)}
.MJXp-scale4 {-webkit-transform: scaleX(.4); -moz-transform: scaleX(.4); -ms-transform: scaleX(.4); -o-transform: scaleX(.4); transform: scaleX(.4)}
.MJXp-scale5 {-webkit-transform: scaleX(.5); -moz-transform: scaleX(.5); -ms-transform: scaleX(.5); -o-transform: scaleX(.5); transform: scaleX(.5)}
.MJXp-scale6 {-webkit-transform: scaleX(.6); -moz-transform: scaleX(.6); -ms-transform: scaleX(.6); -o-transform: scaleX(.6); transform: scaleX(.6)}
.MJXp-scale7 {-webkit-transform: scaleX(.7); -moz-transform: scaleX(.7); -ms-transform: scaleX(.7); -o-transform: scaleX(.7); transform: scaleX(.7)}
.MJXp-scale8 {-webkit-transform: scaleX(.8); -moz-transform: scaleX(.8); -ms-transform: scaleX(.8); -o-transform: scaleX(.8); transform: scaleX(.8)}
.MJXp-scale9 {-webkit-transform: scaleX(.9); -moz-transform: scaleX(.9); -ms-transform: scaleX(.9); -o-transform: scaleX(.9); transform: scaleX(.9)}
.MathJax_PHTML .noError {vertical-align: ; font-size: 90%; text-align: left; color: black; padding: 1px 3px; border: 1px solid}
</style><style type="text/css">.MathJax_Display {text-align: center; margin: 0; position: relative; display: block!important; text-indent: 0; max-width: none; max-height: none; min-width: 0; min-height: 0; width: 100%}
.MathJax .merror {background-color: #FFFF88; color: #CC0000; border: 1px solid #CC0000; padding: 1px 3px; font-style: normal; font-size: 90%}
.MathJax .MJX-monospace {font-family: monospace}
.MathJax .MJX-sans-serif {font-family: sans-serif}
#MathJax_Tooltip {background-color: InfoBackground; color: InfoText; border: 1px solid black; box-shadow: 2px 2px 5px #AAAAAA; -webkit-box-shadow: 2px 2px 5px #AAAAAA; -moz-box-shadow: 2px 2px 5px #AAAAAA; -khtml-box-shadow: 2px 2px 5px #AAAAAA; filter: progid:DXImageTransform.Microsoft.dropshadow(OffX=2, OffY=2, Color='gray', Positive='true'); padding: 3px 4px; z-index: 401; position: absolute; left: 0; top: 0; width: auto; height: auto; display: none}
.MathJax {display: inline; font-style: normal; font-weight: normal; line-height: normal; font-size: 100%; font-size-adjust: none; text-indent: 0; text-align: left; text-transform: none; letter-spacing: normal; word-spacing: normal; word-wrap: normal; white-space: nowrap; float: none; direction: ltr; max-width: none; max-height: none; min-width: 0; min-height: 0; border: 0; padding: 0; margin: 0}
.MathJax:focus, body :focus .MathJax {display: inline-table}
.MathJax.MathJax_FullWidth {text-align: center; display: table-cell!important; width: 10000em!important}
.MathJax img, .MathJax nobr, .MathJax a {border: 0; padding: 0; margin: 0; max-width: none; max-height: none; min-width: 0; min-height: 0; vertical-align: 0; line-height: normal; text-decoration: none}
img.MathJax_strut {border: 0!important; padding: 0!important; margin: 0!important; vertical-align: 0!important}
.MathJax span {display: inline; position: static; border: 0; padding: 0; margin: 0; vertical-align: 0; line-height: normal; text-decoration: none}
.MathJax nobr {white-space: nowrap!important}
.MathJax img {display: inline!important; float: none!important}
.MathJax * {transition: none; -webkit-transition: none; -moz-transition: none; -ms-transition: none; -o-transition: none}
.MathJax_Processing {visibility: hidden; position: fixed; width: 0; height: 0; overflow: hidden}
.MathJax_Processed {display: none!important}
.MathJax_ExBox {display: block!important; overflow: hidden; width: 1px; height: 60ex; min-height: 0; max-height: none}
.MathJax .MathJax_EmBox {display: block!important; overflow: hidden; width: 1px; height: 60em; min-height: 0; max-height: none}
.MathJax_LineBox {display: table!important}
.MathJax_LineBox span {display: table-cell!important; width: 10000em!important; min-width: 0; max-width: none; padding: 0; border: 0; margin: 0}
.MathJax .MathJax_HitBox {cursor: text; background: white; opacity: 0; filter: alpha(opacity=0)}
.MathJax .MathJax_HitBox * {filter: none; opacity: 1; background: transparent}
#MathJax_Tooltip * {filter: none; opacity: 1; background: transparent}
@font-face {font-family: MathJax_Main; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Main-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Main-Regular.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Main-bold; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Main-Bold.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Main-Bold.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Main-italic; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Main-Italic.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Main-Italic.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Math-italic; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Math-Italic.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Math-Italic.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Caligraphic; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Caligraphic-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Caligraphic-Regular.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Size1; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Size1-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Size1-Regular.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Size2; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Size2-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Size2-Regular.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Size3; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Size3-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Size3-Regular.otf?V=2.7.1') format('opentype')}
@font-face {font-family: MathJax_Size4; src: url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/woff/MathJax_Size4-Regular.woff?V=2.7.1') format('woff'), url('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/fonts/HTML-CSS/TeX/otf/MathJax_Size4-Regular.otf?V=2.7.1') format('opentype')}
.MathJax .noError {vertical-align: ; font-size: 90%; text-align: left; color: black; padding: 1px 3px; border: 1px solid}
</style></head>
<body style=""><div style="visibility: hidden; overflow: hidden; position: absolute; top: 0px; height: 1px; width: auto; padding: 0px; border: 0px; margin: 0px; text-align: left; text-indent: 0px; text-transform: none; line-height: normal; letter-spacing: normal; word-spacing: normal;"><div id="MathJax_Hidden"></div></div><div id="MathJax_Message" style="display: none;"></div>
  <div tabindex="-1" id="notebook" class="border-box-sizing">
    <div class="container" id="notebook-container">

<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<h1 id="Formulae-Auto-Drivation-for-Nested-Logit-Models">Formulae Auto-Drivation for Nested Logit Models<a class="anchor-link" href="http://localhost:8888/nbconvert/html/mdm_model_formula_check.ipynb?download=false#Formulae-Auto-Drivation-for-Nested-Logit-Models">¶</a></h1><p>This file uses the Symbolic Python (SymPy) package to derive the elasticity formulae for Nested Logit Models. SymPy is a Python library for symbolic mathematics. It aims to become a full-featured computer algebra system (CAS) while keeping the code as simple as possible in order to be comprehensible and easily extensible. More details about SymPy can be found on the SymPy package main page: <a href="https://www.sympy.org/en/index.html">SymPy</a>.</p>
<p>To use the SymPy package, we first have to import the SymPy package as follows. Note that we also include the <code>init_printing()</code> function to enable the pretty print which is particularly useful when we are printing a large amount of mathematical formulae.</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[1]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#### Import libraries</span>
<span class="kn">from</span> <span class="nn">sympy</span> <span class="k">import</span> <span class="o">*</span>
<span class="kn">import</span> <span class="nn">sys</span>
<span class="n">init_printing</span><span class="p">()</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>Next, let's set up our model by specifying a few basic parameters in the model.</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[2]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#### Settings</span>
<span class="n">number_end_group</span> <span class="o">=</span> <span class="mi">3</span>  <span class="c1"># This is the number of vehicles in each end group of our tree</span>
<span class="n">number_splits_non_end_node</span> <span class="o">=</span> <span class="mi">2</span>  <span class="c1"># This is the number of splits for each node if the node is not on the last two levels</span>
<span class="n">depth_of_tree</span> <span class="o">=</span> <span class="mi">3</span>  <span class="c1"># This is the depth of our tree</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>From the above settings, we build up a tree with the following structre.</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[49]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="kn">import</span> <span class="nn">matplotlib.pyplot</span> <span class="k">as</span> <span class="nn">plt</span>
<span class="o">%</span><span class="k">matplotlib</span> inline

<span class="n">decision_node</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">(</span><span class="n">boxstyle</span><span class="o">=</span><span class="s2">"sawtooth"</span><span class="p">,</span> <span class="n">fc</span><span class="o">=</span><span class="s2">"0.8"</span><span class="p">)</span>
<span class="n">leaf_node</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">(</span><span class="n">boxstyle</span><span class="o">=</span><span class="s2">"round4"</span><span class="p">,</span> <span class="n">fc</span><span class="o">=</span><span class="s2">"0.8"</span><span class="p">)</span>
<span class="n">arrow_args</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">(</span><span class="n">arrowstyle</span><span class="o">=</span><span class="s2">"&lt;-"</span><span class="p">)</span>


<span class="k">def</span> <span class="nf">get_leaf_num</span><span class="p">(</span><span class="n">tree</span><span class="p">):</span>
    <span class="n">leaf_num</span> <span class="o">=</span> <span class="mi">0</span>
    <span class="k">for</span> <span class="n">key</span> <span class="ow">in</span> <span class="n">tree</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span>
        <span class="k">if</span> <span class="nb">type</span><span class="p">(</span><span class="n">tree</span><span class="p">[</span><span class="n">key</span><span class="p">])</span><span class="o">.</span><span class="vm">__name__</span> <span class="o">==</span> <span class="s2">"dict"</span><span class="p">:</span>
            <span class="n">leaf_num</span> <span class="o">+=</span><span class="n">get_leaf_num</span><span class="p">(</span><span class="n">tree</span><span class="p">[</span><span class="n">key</span><span class="p">])</span>
        <span class="k">else</span><span class="p">:</span>
            <span class="n">leaf_num</span> <span class="o">+=</span><span class="mi">1</span>
    <span class="k">return</span> <span class="n">leaf_num</span>


<span class="k">def</span> <span class="nf">get_tree_depth</span><span class="p">(</span><span class="n">tree</span><span class="p">):</span>
    <span class="n">depth</span> <span class="o">=</span> <span class="mi">0</span>
    <span class="k">for</span> <span class="n">key</span> <span class="ow">in</span> <span class="n">tree</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span>
        <span class="k">if</span> <span class="nb">type</span><span class="p">(</span><span class="n">tree</span><span class="p">[</span><span class="n">key</span><span class="p">])</span><span class="o">.</span><span class="vm">__name__</span> <span class="o">==</span> <span class="s2">"dict"</span><span class="p">:</span>
            <span class="n">thisdepth</span> <span class="o">=</span> <span class="mi">1</span><span class="o">+</span> <span class="n">get_tree_depth</span><span class="p">(</span><span class="n">tree</span><span class="p">[</span><span class="n">key</span><span class="p">])</span>
        <span class="k">else</span><span class="p">:</span>
            <span class="n">thisdepth</span> <span class="o">=</span> <span class="mi">1</span>
        <span class="k">if</span> <span class="n">thisdepth</span><span class="o">&gt;</span><span class="n">depth</span><span class="p">:</span> <span class="n">depth</span> <span class="o">=</span> <span class="n">thisdepth</span>
    <span class="k">return</span> <span class="n">depth</span>


<span class="k">def</span> <span class="nf">plotNode</span><span class="p">(</span><span class="n">nodeTxt</span><span class="p">,</span> <span class="n">centerPt</span><span class="p">,</span> <span class="n">parentPt</span><span class="p">,</span> <span class="n">nodeType</span><span class="p">):</span>
    <span class="n">createPlot</span><span class="o">.</span><span class="n">ax1</span><span class="o">.</span><span class="n">annotate</span><span class="p">(</span><span class="n">nodeTxt</span><span class="p">,</span> <span class="n">xy</span><span class="o">=</span><span class="n">parentPt</span><span class="p">,</span> <span class="n">xycoords</span><span class="o">=</span><span class="s1">'axes fraction'</span><span class="p">,</span>
                            <span class="n">xytext</span><span class="o">=</span><span class="n">centerPt</span><span class="p">,</span> <span class="n">textcoords</span><span class="o">=</span><span class="s1">'axes fraction'</span><span class="p">,</span>
                            <span class="n">va</span><span class="o">=</span><span class="s2">"center"</span><span class="p">,</span> <span class="n">ha</span><span class="o">=</span><span class="s2">"center"</span><span class="p">,</span> <span class="n">bbox</span><span class="o">=</span><span class="n">nodeType</span><span class="p">,</span> <span class="n">arrowprops</span><span class="o">=</span><span class="n">arrow_args</span><span class="p">)</span>

    
<span class="k">def</span> <span class="nf">plotTree</span><span class="p">(</span><span class="n">myTree</span><span class="p">,</span> <span class="n">parentPt</span><span class="p">,</span> <span class="n">nodeTxt</span><span class="p">):</span>
    <span class="n">numLeafs</span> <span class="o">=</span> <span class="n">get_leaf_num</span><span class="p">(</span><span class="n">myTree</span><span class="p">)</span>
    <span class="n">depth</span> <span class="o">=</span> <span class="n">get_tree_depth</span><span class="p">(</span><span class="n">myTree</span><span class="p">)</span>
    <span class="n">firstStr</span> <span class="o">=</span> <span class="nb">list</span><span class="p">(</span><span class="n">myTree</span><span class="o">.</span><span class="n">keys</span><span class="p">())[</span><span class="mi">0</span><span class="p">]</span>
    <span class="n">cntrPt</span> <span class="o">=</span> <span class="p">(</span><span class="n">plotTree</span><span class="o">.</span><span class="n">xOff</span> <span class="o">+</span> <span class="p">(</span><span class="mf">1.0</span> <span class="o">+</span> <span class="nb">float</span><span class="p">(</span><span class="n">numLeafs</span><span class="p">))</span> <span class="o">/</span> <span class="mf">2.0</span> <span class="o">/</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">totalW</span><span class="p">,</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span><span class="p">)</span>
    <span class="n">plotNode</span><span class="p">(</span><span class="n">firstStr</span><span class="p">,</span> <span class="n">cntrPt</span><span class="p">,</span> <span class="n">parentPt</span><span class="p">,</span> <span class="n">decision_node</span><span class="p">)</span>
    <span class="n">secondDict</span> <span class="o">=</span> <span class="n">myTree</span><span class="p">[</span><span class="n">firstStr</span><span class="p">]</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span> <span class="o">=</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span> <span class="o">-</span> <span class="mf">1.0</span> <span class="o">/</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">totalD</span>
    <span class="k">for</span> <span class="n">key</span> <span class="ow">in</span> <span class="n">myTree</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span>
        <span class="k">if</span> <span class="nb">type</span><span class="p">(</span><span class="n">myTree</span><span class="p">[</span><span class="n">key</span><span class="p">])</span><span class="o">.</span><span class="vm">__name__</span> <span class="o">==</span> <span class="s1">'dict'</span><span class="p">:</span>
            <span class="n">plotTree</span><span class="p">(</span><span class="n">myTree</span><span class="p">[</span><span class="n">key</span><span class="p">],</span> <span class="n">cntrPt</span><span class="p">,</span> <span class="nb">str</span><span class="p">(</span><span class="n">key</span><span class="p">))</span>
        <span class="k">else</span><span class="p">:</span>
            <span class="n">plotTree</span><span class="o">.</span><span class="n">xOff</span> <span class="o">=</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">xOff</span> <span class="o">+</span> <span class="mf">1.0</span> <span class="o">/</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">totalW</span>
            <span class="n">plotNode</span><span class="p">(</span><span class="n">myTree</span><span class="p">[</span><span class="n">key</span><span class="p">],</span> <span class="p">(</span><span class="n">plotTree</span><span class="o">.</span><span class="n">xOff</span><span class="p">,</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span><span class="p">),</span> <span class="n">cntrPt</span><span class="p">,</span> <span class="n">leaf_node</span><span class="p">)</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span> <span class="o">=</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span> <span class="o">+</span> <span class="mf">1.0</span> <span class="o">/</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">totalD</span>


<span class="k">def</span> <span class="nf">createPlot</span><span class="p">(</span><span class="n">inTree</span><span class="p">):</span>
    <span class="n">fig</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">facecolor</span><span class="o">=</span><span class="s1">'white'</span><span class="p">)</span>
    <span class="n">fig</span><span class="o">.</span><span class="n">clf</span><span class="p">()</span>
    <span class="n">axprops</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">(</span><span class="n">xticks</span><span class="o">=</span><span class="p">[],</span> <span class="n">yticks</span><span class="o">=</span><span class="p">[])</span>
    <span class="n">createPlot</span><span class="o">.</span><span class="n">ax1</span> <span class="o">=</span> <span class="n">plt</span><span class="o">.</span><span class="n">subplot</span><span class="p">(</span><span class="mi">111</span><span class="p">,</span> <span class="n">frameon</span><span class="o">=</span><span class="kc">False</span><span class="p">,</span> <span class="o">**</span><span class="n">axprops</span><span class="p">)</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">totalW</span> <span class="o">=</span> <span class="nb">float</span><span class="p">(</span><span class="n">get_leaf_num</span><span class="p">(</span><span class="n">inTree</span><span class="p">))</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">totalD</span> <span class="o">=</span> <span class="nb">float</span><span class="p">(</span><span class="n">get_tree_depth</span><span class="p">(</span><span class="n">inTree</span><span class="p">))</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">xOff</span> <span class="o">=</span> <span class="o">-</span><span class="mf">0.5</span> <span class="o">/</span> <span class="n">plotTree</span><span class="o">.</span><span class="n">totalW</span>
    <span class="n">plotTree</span><span class="o">.</span><span class="n">yOff</span> <span class="o">=</span> <span class="mf">1.0</span>
    <span class="n">plotTree</span><span class="p">(</span><span class="n">inTree</span><span class="p">,</span> <span class="p">(</span><span class="mf">0.5</span><span class="p">,</span> <span class="mf">1.0</span><span class="p">),</span> <span class="s1">''</span><span class="p">)</span>
    <span class="n">plt</span><span class="o">.</span><span class="n">show</span><span class="p">()</span>


<span class="k">def</span> <span class="nf">create_dict</span><span class="p">(</span><span class="n">level</span><span class="p">):</span>
    <span class="k">global</span> <span class="n">count</span><span class="p">,</span> <span class="n">number_end_group</span><span class="p">,</span> <span class="n">number_splits_non_end_node</span><span class="p">,</span> <span class="n">depth_of_tree</span>
    <span class="k">if</span> <span class="n">level</span> <span class="o">==</span> <span class="n">depth_of_tree</span> <span class="o">-</span> <span class="mi">1</span><span class="p">:</span>
        <span class="n">dic</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">()</span>
        <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">number_end_group</span><span class="p">):</span>
            <span class="n">dic</span><span class="p">[</span><span class="nb">str</span><span class="p">(</span><span class="n">i</span><span class="p">)]</span> <span class="o">=</span> <span class="n">count</span>
            <span class="n">count</span> <span class="o">+=</span> <span class="mi">1</span>
        <span class="k">return</span> <span class="n">dic</span>
    <span class="k">else</span><span class="p">:</span>
        <span class="n">dic</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">()</span>
        <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">number_splits_non_end_node</span><span class="p">):</span>
            <span class="n">name</span> <span class="o">=</span> <span class="nb">str</span><span class="p">(</span><span class="n">i</span><span class="p">)</span>
            <span class="n">dic</span><span class="p">[</span><span class="nb">str</span><span class="p">(</span><span class="n">i</span><span class="p">)]</span> <span class="o">=</span> <span class="n">create_dict</span><span class="p">(</span><span class="n">level</span><span class="o">+</span><span class="mi">1</span><span class="p">)</span>
        <span class="k">return</span> <span class="n">dic</span>

    
<span class="n">count</span> <span class="o">=</span> <span class="mi">0</span>
<span class="n">createPlot</span><span class="p">(</span><span class="n">create_dict</span><span class="p">(</span><span class="mi">0</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt"></div>




<div class="output_png output_subarea ">
<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAWQAAADxCAYAAAD8x81kAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADl0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uIDIuMi4yLCBodHRwOi8vbWF0cGxvdGxpYi5vcmcvhp/UCwAAIABJREFUeJzs3WdcFNf7NvCLKt1GE8ReaAqKDRBEKSooCGIHC4q/mNgribFrbBjsXUAFVFRsIBYsWLChiHSxgyAIiIB09jwvfNh/TEBXtswOe76vEtg953I/y72zM2fuI0UIIaAoiqIYJ810AIqiKOorWpApiqLEhMzKlStXMh2Conh16NAh9OvXD+3bt0daWhp69OiBpk2bori4GMbGxigrK4OKigqMjIyQmZmJoUOHQkpKiunYFMUTKXoOmWKLa9euYdy4cViyZAk2b96Mmpoa+Pj4YNu2bSguLsbvv/+OQ4cOISsrC0uWLMGxY8cwadIkLF68mOnoFMUTWaYDUBSvWrVqBUIIFBQUcPjwYXA4HGhqasLExARlZWXQ0dFBnz59UFhYCADIy8tDly5dGE5NUbyjR8gUq2zduhV79uxBSEjIdx+3Zs0aKCsrIzg4WETJKIp/9KIexRoJCQlYu3Yt5s+f/83PY2Ji4ObmhhEjRiAwMBAAMHXqVFy/fh2hoaEMJKWohqGnLCjWyM3NhYyMDHR0dLg/q6mpwcaNG7Fr1y5oaWlh4sSJsLa2RqtWraCuro43b94wF5iifhI9ZUGxyurVq3HmzBkcPHgQAPDs2TPs378fO3fuBAAEBAQAALKyspCUlITU1FS6yoJiDXrKgmKNvLw8BAcHw8bGhvuz3NxcaGlpcf9fU1MTubm5sLW1RUFBARISEhhISlENQwsyxRpXrlxBSUkJ3N3dv/s4KSkp9O3bF/r6+vD39xdROoriHz1lQbEGh8PB+PHj8fnzZ6xduxYZGRlYvXo15OXlsWvXLmzYsAFVVVVo3bo1pKSkcPPmTdy9exdqampMR6contAjZIo1CCGorKyErKwsnj17hmnTpsHOzg6ZmZl4//49BgwYgPDwcACAnJwcCgoKUFNTw3BqiuIdPUKmWCMoKAhLlizBr7/+Cl9fX6xcuRL9+/fHnTt38Pfff6Ompgb9+/fH7du3MWzYMKSkpKBHjx7Ytm0b09Epiie0IFOsUVBQAENDQ3z58gV79+6Fvr5+nY/Ly8uDt7c3cnNzcf/+fZiamoo4KUU1DD1lQbFCdXU1li5dCkII2rVrV28xBgB1dXVYWFigWbNmWLBgAfdWaooSd7QgU2KvuLgYzs7OePToETgcDjZs2IDbt2/j5s2bAIDY2FhERkYCABITE3HmzBnMmjULWlpaKC8vh6WlJb1BhGIFeqceJdbev38PJycn9OnTBxs3boSdnR38/f0RExMDaWlp3L9/Hzdu3ICqqiqePHmC6OhoaGho4OnTp8jKykJgYCBevnwJS0tLnDt3Dr169WL6n0RR9aJHyJTYio+Ph7m5OcaNG4d9+/ahW7duuHr1KioqKnDz5k3cunULZWVluHTpEmJiYsDhcBAaGop79+5BSUkJBw4cgJOTE2bPno3du3fD0dER586dY/qfRVH1ohf1KLF06dIleHp6YufOnRgzZoxAxnz06BFGjBiBxYsXY86cOQIZk6IEiRZkSuzs378fy5cvx+nTp2FpaSnQsd+8eQNHR0fY29vj77//hoyMjEDHpyh+0IJMiQ0Oh4M//vgDp0+fxsWLF9G5c2ehzFNYWIiRI0dCVVUVwcHBUFZWFso8FPWz6DlkSiyUl5dj3LhxuH37Nu7duye0YgwAzZo1Q2RkJJo1awYbGxt8+PBBaHNR1M+gBZliXF5eHuzs7CAlJYVr165BXV1d6HPKy8sjICAAzs7OMDc3R3JystDnpKgfoQWZYlR6ejrMzc1hZWWFkJAQKCgoiGxuKSkpLFu2DKtXr4aNjQ2uX78usrkpqi60IFOMuXv3LqysrLB48WKsX78e0tLMvB09PT0RGhqKcePG4fDhw4xkoCiAXtSjGHLixAnMmjULR48exeDBg5mOAwBISUmBk5MTPD09sXLlSrrTCCVytCBTIkUIwaZNm7Br1y5cuHABJiYmTEf6Rk5ODpydndG1a1ccPHgQ8vLyTEeiJAgtyJTIVFVV4bfffsOjR48QHh4OXV1dpiPVqbS0FB4eHvj06RPCwsLQvHlzpiNREoKeQ6ZEoqioCMOHD0dmZiZu3boltsUYAJSUlHDy5En07NkTFhYWeP36NdORKAlBCzIldJmZmejfvz/at2+P8+fPQ1VVlelIPyQjI4MtW7Zg5syZsLS0xMOHD5mOREkAWpApoXr69CnMzc3h6emJ3bt3Q1aWXQ0Gf/vtN+zfvx/Dhg3DmTNnmI5DNXL0HDIlNJGRkZg0aRJ27979w52ixd3jx4/h4uKCBQsWYO7cuXQFBiUUtCBTQrF3716sWrUKYWFhMDc3ZzqOQLx79w6Ojo4YOHAgtm7dShsTUQJHCzIlUBwOBz4+Pjh37hwuXryIjh07Mh1JoD5//gx3d3coKCjg2LFjUFFRYToS1YjQc8iUwJSVlWHMmDG4f/8+YmJiGl0xBoCmTZvi4sWL0NTUxIABA5Cdnc10JKoRoQWZEoiPHz/C1tYW8vLyuHr1Klq2bMl0JKGRk5PDwYMH4ebmBnNzcyQmJjIdiWokaEGm+JaWlgZzc3PY2toiKCgITZo0YTqS0ElJSWHp0qX466+/MGjQIERFRTEdiWoE6Dlkii+3b9+Gu7s71q9fDy8vL6bjMOLWrVsYNWqURL8GlGDQgkw1WEhICObOnYuQkBDY2dkxHYdRaWlpcHJywrhx47B69Wq6LI5qEFqQqZ9GCMH69euxb98+REREwNjYmOlIYuHjx49wdnZGhw4d4O/vLxGnbijBoueQqZ9SVVUFb29vnD59Gvfu3aPF+B80NDRw/fp1VFRUwMHBAQUFBUxHoliGFmSKZ58/f4aTkxNycnIQHR0NHR0dpiOJHUVFRYSGhqJv376wsLDAy5cvmY5EsQgtyBRPMjIy0L9/f3Tp0gVnz56lN0R8h7S0NDZt2oQ5c+agf//+uH//PtORKJagBZn6oSdPnsDc3BxTpkzBjh076C3DPJoxYwYOHToEZ2dnnD59muk4FAvQi3rUd4WHh8PLywv79u2Dq6sr03FYKS4uDs7Ozpg7dy7mz59PV2BQ9aIFmarX7t27sXbtWpw5cwZ9+/ZlOg6rZWZmwsnJCZaWlti+fTvr2pBSokELMvUfHA4HixYtwsWLFxEREYEOHTowHalRKCoqwujRoyEjI4MTJ07Q8/DUf9BzyNQ3SktLMWrUKDx+/BgxMTG0GAuQmpoaLly4AB0dHVhbWyMrK4vpSJSYoQWZ4srNzcWgQYOgpKSEy5cv0809hUBOTg779+/H6NGjYW5ujoSEBKYjUWKEFmQKAJCamop+/fph8ODBOHLkCL3LTIikpKTg4+ODTZs2wdbWFleuXGE6EiUm6DlkCtHR0Rg9ejQ2btyIyZMnMx1Hoty5cwfu7u5Yu3Ytpk2bxnQcimG0IEu4oKAgLFiwACEhIbC1tWU6jkRKT0+Ho6MjRo0ahbVr10Jamn5xlVS0IEsoQgjWrl2LQ4cOISIiAkZGRkxHkmh5eXlwcXFBmzZtEBAQAAUFBaYjUQygH8USqLKyEl5eXjh//jzu379Pi7EYUFdXx7Vr18DhcGBvb4/8/HymI1EMoAVZwhQWFmLo0KEoKCjAzZs3oa2tzXQk6v+r3TjV0tIS5ubmePHiBdORKBGjBVmCvH37FpaWljA2NkZYWBiUlZWZjkT9i7S0NDZs2ICFCxfCysoKMTExTEeiRIgWZAkRGxsLCwsLTJ8+Hdu2baMNgsTc9OnTERAQgBEjRuDkyZNMx6FEhF7UkwDnz5/H1KlTcfDgQbi4uDAdh/oJ8fHxGDZsGGbNmoVFixbRxkSNHC3IjdyOHTuwfv16nDt3Dr1792Y6DtUAmZmZGDZsGPr164edO3fSxkSNGC3IjVRNTQ0WLlyIy5cv4+LFi2jXrh3TkSg+FBcXY/To0QCA0NBQqKqqMpyIEgZ6DrkRKi0thbu7O+Lj43H37l1ajBsBVVVVXLhwAW3btoWVlRUyMzOZjkQJAS3IjUxOTg5sbGygpqaGS5cu0QZBjYisrCz27NmDCRMmwMLCAvHx8UxHogSMFuRGJDk5Gf369cOwYcMQGBgIeXl5piNRAiYlJYVFixZhy5YtsLe3R2RkJNORKAGi55AbiRs3bmDs2LHYvHkzJk6cyHQcSgRiYmIwcuRIrFq1CtOnT2c6DiUAtCA3AkeOHMGiRYtw/PhxDBw4kOk4lAi9ePECjo6OcHV1xfr162ljIpajBZnFCCFYtWoVjhw5goiICBgYGDAdiWJAfn4+RowYgVatWuHw4cNQVFRkOhLVQPTjlKUqKysxadIkXLx4Effu3aPFWIK1bNkSV69ehYyMDOzs7JCXl8d0JKqBaEFmoU+fPmHIkCEoLi7GzZs3oaWlxXQkimEKCgoIDg6GjY0NzM3NkZ6eznQkqgFoQWaZ169fw8LCAiYmJjh16hSUlJSYjkSJCWlpaaxbtw5LliyBlZUV7ty5w3Qk6ifRgswiDx8+hKWlJX799Vf4+fnRBkFUnaZNm4YjR47Azc0Nx48fZzoO9RNoQRZj0dHR3K+eZ8+ehZOTE/bu3YtZs2YxnIwSdw4ODoiKisKSJUuwfv16EEJQVFSE0NBQpqNR30FXWYgpDocDfX19BAQE4NGjR9i8eTPOnz8PMzMzpqNRLJKVlQUnJyf06tUL69atQ9euXZGeng51dXWmo1F1oAVZTF29ehULFy6EtbU1bty4gYiICLRt25bpWBQLFRcXY+zYsaiuroampiaMjY2xZMkSpmNRdaCnLISMw+Fg+/btiI2NBSEEe/fuRXR0NAAgODgY4eHhAL6ekvhnI/Lt27eDEIKEhASsXLkSFy5cYCQ/xX6RkZHw9PSEjo4O7t27h127dqGmpgbA17v9du3aBUIInj59ii1btqCmpgbp6elYt24dKioq8P79e6xevRpfvnxh+F/S+NEjZCHicDj49ddfcePGDeTn52PIkCG4f/8+CgsLMXz4cERFRaGiogLOzs4IDw+HtLQ0VqxYASsrKxgZGaF58+aorKyEmZkZhg0bhoULF9IG5dRPO3HiBI4ePYrbt29DWVkZ2dnZ2LNnDwwNDeHm5obmzZujV69eiIqKgoaGBjp37oxHjx6hZcuW0NbWxsuXL6GkpIQWLVogMjKSbv0lRLQgC9GnT5+gqamJffv2obS0FJGRkZg/fz4yMjIQEhKCefPmobi4GPv378evv/6Ke/fu4fz58zh16hT+/PNPzJw5E5aWlnRpGyUQ1dXVePz4Mfbs2QN7e3uEh4fjy5cvmD9/Pnx9fWFnZ4fevXvDz88PpqamsLe3x/bt26Gnpwc7Ozu4uroiIiIC1tbWTP9TGi1akIUsJCQECxYsQEhICNTU1Op9XGpqKmbPno2rV6+iZ8+eIkxISarc3FxYWFhg4sSJcHJy+u5jp0+fjoEDB8LX11dE6SQTPYcsZJ8/f4acnNwPm77UbstDz9NRolJeXo6KigooKCj88LGKioooLCwEPX4TLnqELET5+fnQ1NTE4cOHv+k1ERMTA19fX3A4HIwYMQKTJ08GAISFheHYsWN49eoVQ4kpSeLu7g5CCP744w/uz+p7b5aUlGD48OE4c+YMbG1tGUrc+NEjZCFq0aIFZsyYAT8/P1RUVAD4utfdxo0bsX37dpw8eRKXL1/Gq1evkJOTg8DAQCxfvpzh1JSkWLRoEW7evImnT58CqP+9CQA7d+6EmZkZLC0tmYzc6NGCLERSUlIYP348UlNTUVJSAgBISkqCnp4eWrduDTk5OTg4OCA6Ohrv3r1DVVUVHBwcGE5NSQoTExN06NABCQkJAOp/bwJAbGws3N3deTq9QTUcLchCVFRUhGHDhmHdunVo2bIlgK8XUv7ZnU1TUxO5ubno3bs33NzcMGLECKbiUhJm/vz5UFVVxfjx4wHU/94EAF9fX/z555948uQJI1klBS3IQqSgoAADAwM8fPjwuxdDpKSkUFFRgbi4OHprNCUyZmZmSE9PR0FBQb2PqV33npiYCBUVFWhra4sqnkSiBVmI5OXl4e/vj+PHjyMnJwfA1/N0t2/fRnV1NaZNm4bMzExoaGjg7t27eP36Nfz8/BhOTUmKqVOnonPnztw7RDU0NHD79m28fPkSu3fvxv3796GhoQHg652jvr6+0NHRYTJyo0cLshCVl5dj9OjR8Pb2hra2NgghuHDhAmpqapCTkwNNTU2cPn0a1tbWsLa2hpGREaZNm8Z0bEpCbN26FRkZGZgwYQIA4OPHjygqKoKcnBzMzMxw9epV9OvXDwDw+++/Y/bs2Xj79i2TkRs9WaYDNGbV1dXIycnhNgW6e/cusrOzsWzZMsyaNQuVlZX48uULFBUVISsrCz09PfqGp0QmIyMDmpqaUFJSQkVFBXbs2IHp06dj7ty5qKmpQdu2bXH37l0YGBigVatWqKqqQlFREdOxGzW6DlnI4uLiYG5ujpMnT2LmzJlYsGDBN0uHDh48iPT0dLi4uGDNmjVISUnhXgCkKGGqrq6GtbU1evToASkpKSQlJX1zJ15WVhY8PT0REhICb29vrFmzBl5eXgwmbvzoKQsh4nA4+Pvvv9G7d29ERUVBV1eXW4yrq6sBAB4eHkhKSkJlZSWUlZURFBTEZGRKgkRHR+P58+cwMDBAUFAQ5syZA+DrdQ5CCHR0dODm5oYdO3bAwcEBu3fvRnFxMcOpGzdakIXo8+fPOHHiBJydnbFv3z7ExcUhNjYWBw8ehJWVFW7evInw8HDk5ORg8+bNcHZ2xu7du5mOTUmIQ4cOoU+fPjh27BhKS0sRFBSE9PR0DB06FEuXLsXbt29x6dIlXL9+HV27dkVaWhr3JhJKSAglVOfPnycyMjJES0uLhIWFETU1NaKvr0/Cw8OJuro60dPTIxEREUReXp6oqqqSxMREpiNTEiI/P5907tyZSElJkSNHjhALCwvSpEkTEhAQQBwdHYmcnBzZtm0bsbCwIADI6tWrmY7c6NGCLGRPnjwh6urqJCMjgxBCyJs3b8inT58IIYS8f/+efPz4kRBCSFRUFNHQ0OD+jqKEjcPhkN69e5PNmzcTQgj58uULef78OSGEkIqKCpKcnEwIIaS6upqYmJiQgIAApqJKDHpRT4gIIRgwYAA8PDwwffr0Hz5++vTpUFFRwd9//y2CdJSkCw4OxtatW/HgwYMfdiN89OgRXFxckJqa+t02shR/aEEWotDQUPz11194/PgxZGRkfvj43NxcGBoa4s6dO9DX1xdBQkpSlZSUQF9fH6GhobCwsODpOZMnT4a2tjY2bNgg5HSSixZkISkrK4OBgQECAwNhY2PD8/O2bNmC69evIyIiQnjhKIm3bNkyvHr1CsHBwTw/Jzs7G926dcODBw/QsWNHIaaTXLQgC8maNWvw7NmzbzYu5UVlZSWMjY2xdetWODo6CikdJcnevHkDMzMzxMfHo3Xr1j/13A0bNuDBgwc4c+aMkNJJNlqQhSAjIwOmpqZ4/Pgx2rVr99PPj4iIwIIFC/Ds2TPIy8sLPiAl0UaNGoXu3btj2bJlP/3c8vJyGBkZYd++fbCzsxNCOslG1yELgY+PD3799dcGFWMAcHR0RPv27bFr1y7BBqMkXnR0NB49eoSFCxc26PkKCgrw9fXF3LlzuTc3UYJDj5AFLCYmBmPGjEFqaipf26WnpKTA2toaycnJ3I5bFMWPmpoamJmZYenSpRg1alSDxyGEwM7ODm5ubvjtt98EmJCiBVmAOBwO+vbtizlz5sDDw4Pv8ebOnYvy8nLs3btXAOkoSbd//34EBQUhOjqa2+e4oRISEmBnZ4eUlBS0aNFCQAkpWpAFKDAwEPv27cPdu3d/uK6TF58+fYK+vj4uX74MU1NTASSkJFVhYSH09fURGRmJHj16CGTM3377DTIyMti+fbtAxqNoQRaYoqIi6Ovr4+zZs+jTp4/Axt27dy+OHz+OGzdu8H1UQ0mu+fPno6SkBPv37xfYmHl5eTA0NMSNGzdgZGQksHElGS3IAuLj44MPHz4gMDBQoOPW1NSgZ8+eWLZsGdzd3QU6NiUZ0tLS0L9/fyQlJUFTU1OgY2/fvh3h4eG4fPkyPWAQAFqQBeDFixfo168fEhIS0KpVK4GPf+PGDUyZMgUpKSlQVFQU+PhU4+bk5ARbW1vMnz9f4GNXVVXBxMQEGzZsgLOzs8DHlzR02ZsALFy4EAsWLBBKMQaAgQMHwszMDFu2bBHK+FTjdfHiRbx48QIzZ84UyvhycnLw8/PD/PnzUVFRIZQ5JAk9QuZTVFQUpk+fjuTkZCgoKAhtnlevXqF379549uwZdHV1hTYP1XhUVlaie/fu2LJlC5ycnIQ6l7OzM6ysrLBo0SKhztPY0YLMh+rqapiammLNmjVwdXUV+nxLly7Fu3fvcPToUaHPRbGfn58frly5gosXLwr9/G56ejrMzc2RmJgIbW1toc7VmNGCzIddu3YhLCwMUVFRIrmgUVJSgq5du+LUqVMwNzcX+nwUe338+BGGhoa4desWDAwMRDLnokWLUFBQgEOHDolkvsaIFuQGys/Ph4GBAa5du4Zu3bqJbN6jR49ix44duH//vkDWOlON0y+//AJFRUX4+fmJbM7Pnz9DX18fFy5cQK9evUQ2b2NCC3IDzZo1CxwOR+T9JjgcDiwsLDBjxgxMmjRJpHNT7PD06VMMGTIEKSkpaN68uUjnPnToEPz9/XHnzh26DK4BaEFugKSkJNjY2CAlJQXq6uoin//Bgwdwc3NDamoqVFVVRT4/Jb4IIRg4cCDGjh2LX375ReTz19TUoE+fPli4cCHGjRsn8vnZjn7n/UmEEMybNw/Lli1jpBgDQN++fWFra4v169czMj8lvsLCwvDp0yd4e3szMr+MjAy2bt2KJUuWoLS0lJEMbEaPkH/S+fPn4ePjg/j4eMjJyTGW4/379+jevTsePXqEDh06MJaDEh9lZWUwNDSEv78/Bg4cyGiWsWPHQl9fHytXrmQ0B9vQgvwTKioqYGRkhF27dmHw4MFMx8Fff/2F2NhYhIWFMR2FEgPr1q1DXFwcTp06xXQUvHv3Dj169EBcXBzatGnDdBzWoAX5J2zatAm3b9/GhQsXmI4C4OvuDYaGhjhw4ABsbW2ZjkMx6P379zAxMcGjR4/Qvn17puMAAFauXImUlBScOHGC6SisQQsyjz58+ABjY2Pcu3cPnTt3ZjoOV1hYGFasWIG4uDjIysoyHYdiiKenJ9q2bYu1a9cyHYWrtLQU+vr6CA4OhpWVFdNxWIFe1OPR0qVLMWXKFLEqxgDg6uoKdXV1HDhwgOkoFEPu37+PGzduwMfHh+ko31BSUsLGjRsxZ84c1NTUMB2HFegRMg9iY2MxfPhwpKamomnTpkzH+Y/4+Hg4ODjQ3RskEIfDgbm5OWbOnAlPT0+m4/wHIQRWVlaYMmUKpk6dynQcsUcL8g8QQtC/f394eXmJ9RtqxowZkJeXx7Zt25iOQonQkSNHsHv3bsTExIjtnZuPHz/GsGHDxPaARpzQgvwDx44dg6+vLx4+fAgZGRmm49SrtndBdHQ0DA0NmY5DiUBxcTH09fURFhaGvn37Mh3nu6ZOnYoWLVpg8+bNTEcRa7Qgf8eXL1+gr6+PY8eOoX///kzH+aFt27YhIiKC7t4gIf744w+8f/8ehw8fZjrKD9VeFI+JiUGXLl2YjiO2aEH+jhUrVuD58+c4duwY01F4UlVVhe7du2PTpk0YPnw403EoIXr16hX69OmDZ8+eQUdHh+k4PNm8eTNu3bolNstGxREtyPVg68L2S5cuYdasWUhKSoK8vDzTcSghGTlyJMzMzPDHH38wHYVnFRUVMDY2xs6dO8XixipxJJ5XAcTA4sWLMWvWLFYVYwAYMmQIunbtSrdmb8SuX7+OuLg4oeyRJ0xNmjTBli1bMG/ePFRVVTEdRyzRI+Q63L59GxMmTEBqaiqUlJSYjvPTnj9/DgsLCyQlJUFLS4vpOJQAVVdXo2fPnli5ciXc3NyYjvPTCCEYPHgwhg0bhtmzZzMdR+zQgvwvNTU16N27NxYvXoyxY8cyHafBFi5ciMLCQhw8eJDpKJQA7dmzBydPnsS1a9dYe+GW6fa14owW5H85ePAgAgMDcfv2bda+4YGvuzd07doVERERMDMzYzoOJQAFBQUwMDDA1atX0b17d6bj8GX27Nmorq7G7t27mY4iVmhB/ofGVsQOHDiAw4cPs/7Dhfpqzpw5qKysxJ49e5iOwreCggLo6+sjKiqK9R8ugkQL8j80tk0aa2pq0KtXL/j4+GDMmDFMx6H4kJycjAEDBiA5ORkaGhpMxxGIXbt24fTp06w+/SJotCD/f7UXwhrbNua3bt2Cp6cnUlJSWHmBkvp6IWzIkCFwdHTEnDlzmI4jMNXV1TA1NcWaNWvg6urKdByxQJe9/X8LFizAkiVLGlUxBgBra2v069eP3rLKYhEREXj37h1+/fVXpqMIlKysLLZt24YFCxagvLyc6ThigR4h4/9upkhMTESTJk2YjiNwb9++Rc+ePVl3kwsFVFZWwsjIqFHfTOHq6oo+ffrg999/ZzoK4yS+IEvK7cZsuw2c+srX1xfR0dGN+nbjly9fom/fvqy6DVxYJL4gS0pDHrY1SqKAnJwcGBkZSURDHh8fH2RnZ7OiUZIwSXRBzsvLg4GBAW7evAkjIyOm4whdSEgItmzZgkePHolt71zq/3h7e6Np06bw9fVlOorQFRcXo2vXrjh79iz69OnDdBzGSHRB/vXXXyErKysxfR9qm+1PnToVXl5eTMehvuPJkydwcnKSqKbugYGB2Lt3r1g32xc2iS3Iz549g729vcRte1S7HVVaWhrU1NSYjkPVgRACa2trTJo0CdOmTWM6jshwOBz07dsXc+bMgYeHB9NxGCGRH0OEEMydOxcrVqyQqGIMAL169cLQoUPRfctQAAAgAElEQVTFandi6luhoaH48uULpkyZwnQUkZKWlsb27dvh4+ODkpISpuMwQiKPkMPCwrBixQrExcVBVlaW6TgiV7t7w71798RuF21JV1paCgMDAwQFBcHKyorpOIzw8PBAu3btJPKgQWKOkF+/fo1Vq1ahvLwcCxcuxNatWyWyGAOAtrY2Fi1ahIULFwL4evGI9qdl1rx58/Dp0yf4+vqiX79+EluMAWDDhg3Ys2cP3rx5g2PHjuHKlStMRxIZiTlCDg8Px+7du2FlZYWHDx/izJkzTEdiVEVFBYyMjLB79254eHggISGB9k5mkK6uLs6cOQNHR0c8fvwYbdu2ZToSo9asWYNnz55BX18fcnJyWL58OdORREJijpA/f/7M3bFg1qxZGDt2LD5+/Mh0LEbcvXsX8+fPx6pVqzBv3jw0bdoUhYWFTMeSaIWFhdi8eTN++eUXBAUFITAwkOlIjCCEwMPDA5aWlnj06BHy8vLw+fNnpmOJjEQV5KSkJJiammL06NHo3bu3xDbHNjU1BQAsWbIEcnJyqKyslKg3vbipqqpCeXk5bt++jQsXLuDu3buwt7dnOhYjpKSk4O7uDg8PDxgZGeH8+fMSdbAgMQU5MTER6enpKC0tRUxMDBYsWNCo78z7HmVlZezatQtBQUHIz8/Hu3fvkJGRwXQsiVVYWAgOh4OysjLMnz8fERER0NXVZToWY0aMGIGEhAQ0b94c2dnZiIuLYzqSyEjMVa0uXbpg5syZ2LZtm8QuOv83GxsbpKamYuzYsfT8MYMUFRVhYWGB0NBQiS7E/9SyZUsEBQXBwsIC8fHxTMcRGYm5qEdRFCXuWH2oWFVVhZqaGgBfm11XV1cD+LpTRu0yLg6Hg8rKSsYyskl5eTlqP58rKiq4/11ZWQkOhwPg29ecqh8hBBUVFdz//2e/3/r+m6ofL3/f/37N2Yi1BTkzMxNdu3aFpaUlXrx4gR49esDExAQvXryAtbU1OnXqhNTUVAwbNgx6enpISUlhOrJYCw8PR9OmTeHj44Po6Gi0aNEC06ZNQ2xsLLS1tTFy5EgkJSWhXbt2cHBwQGlpKdORxRYhBHPmzEHz5s1x9epVrFixAqqqqggNDcWOHTugrKyMvXv3IiQkBCoqKli/fj3TkcXap0+f0Lt3bxgbGyM9PR02Njbo0KEDUlJS4OLiAh0dHcTFxWHy5MlQV1fHvXv3mI7ccISlDAwMyIwZM8iYMWOItLQ0mTp1Kpk+fTqRlpYm7u7uZNGiRURWVpY4ODiQP/74g2hqapLKykqmY4ullJQUoqqqSrZu3UoMDAyIqqoq2bJlC+nVqxdRUFAga9euJTY2NkROTo78/vvvxN7ennh4eDAdW2xt27aNdOnShfj5+RE1NTXSuXNnsnPnTqKurk5at25Ndu3aRXR1dYmmpibZuXMn0dXVJaGhoUzHFluDBg0iY8aMITNmzCDS0tLE1dWVLFmyhMjKyhJbW1uyevVq0qRJE9K3b1+yYcMGoqamRnJzc5mO3SCsvajXp08fPHv2DL6+vnB3d0f79u0BAA4ODmjbti2kpKTQr18/6OrqYu3atejRowdkZGQYTi2eNDQ0oKWlhYyMDOzduxclJSXQ0tJC3759UVBQAB0dHdja2uLDhw9QVlbGsWPH4OzszHRssWViYoL8/HxoaGjg+PHjUFJSgoqKCoKDgyErK4umTZvi8OHD4HA4SE9PR3l5OfT19ZmOLbb69OmDiIgI7Ny5E7a2tmjTpg2kpaXRt29f6OjoQFZWFiYmJlBXV8e+ffvQpUsXqKioMB27QVh7Ua+oqAh6enr466+/0K9fv3of9/LlS0yYMAGvXr2Cnp6eCBOyS0REBNzc3BATE/PdxwUGBuL+/ft4+PChxC4b5MX//vc/vH79+oenI7y9vTFixAiJuROtIaqqqtCpUyd4eXlh2LBh9T6usLAQdnZ2ePLkCXr06CHChILD2nPI48aNg7W19TfNrGNiYuDm5oYRI0Zw73Tq2LEjJkyYADc3N+5FAepb2dnZ8Pb2xooVK775eV2v56hRo1BRUYGVK1eKPihLnD17FufPn8dvv/3G/VldryXwdXPd7du3//CDUJLNmDEDbdq0+WZPwbpez2bNmmHBggUYPXo0iouLGUrLH9YW5NzcXLRq1Yq7primpgYbN27E9u3bcfLkSVy+fBmvXr0CALRq1Qp5eXlg6ZcBoSsrK0N5eTk0NDS4P6vv9WzSpAlatGiB3NxcBhOLt48fP0JZWZn7tfl7781mzZpBTk4Onz59YjKyWMvJyYGGhga3Gdj3Xk8tLS0UFRWxtlkWawtyeHg4Tpw4wb2LJykpCXp6emjdujXk5OTg4OCA6OhovHv3Dn5+frhy5Qrk5OQYTi2eOnTogAMHDmDBggXcn9X3eoaGhqKoqAjbtm1jMLF48/b2hqWlJfz8/ADU/1oCX5voTJ48GU5OTkxGFmshISF4/PgxoqKiANT/ehYVFWHZsmU4c+YMa/ucs7Yg7969G7q6utx+vrm5ud/cbaapqYnc3Fxoa2vD1NQUmzZtokfI9SgtLcWWLVu+2XW7vtfT0tISmZmZiIiIYCIqKyQmJuLixYtwcHAAUP9rCQBDhw5FUFAQ3r17x0hWNggODgYhhHteuL7XU1lZGTY2Nti4cSNrT0+ytiAfPHgQY8aM+e7VVCkpKcjLy2PixIkIDg5m/aJxYUlNTcXTp08xYcKE7z5OSkoKbdu2ha2tLfz9/UWUjn1Onz6NNm3afPdic+0FUXt7e8jLy+Py5cuiisc6Bw4cgIuLy3ebgUlJSUFGRgaTJ09GZGQkcnJyRJhQcFhbkMPCwrB9+3buuSNNTU28ffsWWVlZyMjIQG5uLjQ0NJCXl4cVK1bg1KlTUFBQYDi1eOrZsyf+/PNPzJ07l/szTU1NvH79Gm/fvkVhYSH39YyKikJMTAz27dvHYGLxtmTJEqiqqnJfI01NTWRnZ+PFixd4/vw595woAGzcuBE9e/akm85+R0hIyDenJzU1NZGRkYHMzEy8efOG+96sbc60f/9+1vYEYe065Pz8fMjIyEBeXh4AYGBggMTERIwdOxZNmjRBy5YtsW7dOsjKykJeXp61n5iiQAhBVlbWN5ueqqmpISEhAXPmzEG3bt3w4sULrF27Fh8/fkR5eTlrr2KLQnl5OQoKCrhtTg0NDZGWlgYfHx/u72vPwauqquL9+/eoqqqi6+Tr8enTJ9TU1EBRURHA19czKSkJ06ZNQ2lpKbS1tbF+/XpISUlBRUUF2dnZDCduONYeIU+aNAlz585F69atAQBpaWlo1qwZysrKUFJSgk6dOqFjx45o1qwZVq5cCW9vb9o3oB4PHz7EoUOHsGHDBu7PLl68iO7duyM3NxeXLl2CtbU1OnbsiH79+mHw4MGYPXs2g4nF2+bNm6GgoMDdOVlGRgbKysrIycnBhw8fIC8vj44dOwIAZs+ejYyMDBw6dIjJyGLN29sbkyZN4t488/HjR8jKyqKwsBA1NTXQ1dVFx44doaCggPXr1+PPP/9k7Tl51hbkNWvWYM+ePcjPzwcAnDt3DqNHj4aioiLmzZvH7bXw5csX+Pn5YenSpfSURT3MzMxgb28PX19fAF+XFYWHh8PZ2RnGxsawsbGBtrY2ACA5ORkRERFYtmwZk5HF2v/+9z9kZ2cjMjISwNeDBQBwdXXF1KlTUVhYyP3GFhISAnl5eYwbN46xvOJu9erVCAoK4vbsvnDhAhwdHaGjo4NFixZxa0B1dTU2b96MqVOnsvYmMNYWZHt7e5SUlODDhw8oKyvDtWvXuHfxDBw4EHFxcdztX969e0dv9f0OWVlZuLq64unTpwCA+/fvQ11dnXsezsXFBefOnQPwtbhoa2vDxMSEsbzirnXr1txb+4GvBwvDhw+HlJQU5OTkYGdnh/DwcABAfHw8Bg0axNplWqJgYWEBeXl5vH37FhwOBxcuXICLiwsAoHv37sjPz+fegp6cnIyRI0ey9i5S1hZkJycneHl5QUtLC3/88QcMDQ2hqKiIqqoqlJaWYsCAAVi+fDlqamqwdOlSDBkyhLWLxYUtNTUVs2fPxrZt2xAbG4tdu3bB2dkZb968QU5ODnr37o33799j48aNcHR0RPv27fHLL78wHVtsbdu2DampqZg7dy6OHTuGixcvYvjw4Xj16hVevHgBFxcXnDhxAhEREfjzzz9x7tw5hIaGMh1bbI0ZMwa2trYwNTXF8uXLoaSkhPbt2+Pz58/IycmBo6MjNmzYgMzMTG5vG9bul8lkZyN+TJ8+nZiampJ27doRNTU1oqurSwwNDYm0tDTR09Mj+vr6REFBgWhpaRELCwvi4uJCampqmI4tlvLz84mhoSEZPHgwadasGZGWliY2NjZETU2NqKiokAEDBhAlJSWira1N+vbtS3R0dMi+ffuYji227ty5Q1q0aEEcHByIjo4OkZOTI05OTkRNTY2oqqoSJycnIisrS7S0tIijoyPR0NAgiYmJTMcWWytXriSdOnUi3bp1I+rq6qRZs2akT58+REFBgbRs2ZL06tWLyMrKkpYtW5JBgwYRS0tLUlZWxnTsBpFZydKmBE5OTsjMzISJiQni4+MxY8YMdOvWDQ8fPoSPjw90dXWRm5sLT09PtGjRAgcOHODeekl9S1FREaNHj0Z0dDQsLS3Rpk0bqKurY+TIkSguLoaBgQHGjRuHW7duYejQoRg7diy8vb2Zji22atcgP3/+HNLS0pg+fTry8/PRq1cvWFlZoaSkBM7OzlBWVoaCggIOHDiAbt26MR1bbNnY2KCiogK6urp4+PAh/vzzT7Rs2RIfPnzA3LlzoaioCCUlJbi4uEBOTg7BwcFQUlJiOnbDMP2JwK/ff/+dzJ8/n/v/KioqpKioiBBCyKZNm8iUKVOYisY6HA6HGBsbkxs3bhBCCLlx4wYZMGAA9/fW1tbk9OnTzIRjoVevXhF1dXVSXl5OCCFk/vz5xNfXlxBCyMePH0nTpk1JYWEhkxFZZceOHWTs2LHc/+/cuTNJS0sjhBASFBREBg8ezFQ0gWHtOWTg61XVw4cP17uofuLEiQgLC6NrZnkUGxuL0tJSWFtb1/n7qVOn0uVZPyEgIADjx49HkyZN/vM7dXV12Nvb4/jx4wwkY6dDhw5h6tSpdf7Ozc0Njx49Yu1yt1qsLsiXL1+Gnp4ejIyM6vy9lpYWbGxs6AUTHh06dAheXl717so9cuRI3Lt3D+/fvxdxMvapqalBYGDgd+/A8/Lyoh9wPIqLi8OnT58waNCgOn+vqKiIMWPG4PDhwyJOJlisLsj+/v71fmLWokd1vCktLUVoaCgmTZpU72OUlZUxatQo1r/pRSEqKgqamprfXR7o4OCArKwsJCQkiDAZOx06dAhTpkyp92AB+Pq37u/vz92Ql41YW5Bzc3Nx7do1jBkz5ruPGzp0KN68eUM3Of2B06dPo1+/ftw7H+vj5eUFf39/2jnvB3g5WKhthkMbNX1feXk5jh8/jsmTJ3/3cT179oSamhpu3rwpklzCwNqCHBQUhBEjRnzTf6EusrKymDhxIn3T/wAvBQT4ur+ZgoICbt26JYJU7JSfn4/Lly/zdPfdlClTEBwczN3KnvqvM2fOoGfPnmjbtu13HyclJcX6b8SsLMiEEO75Tl54eXnhyJEj9MaQerx8+RJJSUnf9EOuj5SUFPcomapbUFAQhg0bhmbNmv3wsR07doSRkRHOnz8vgmTsxOvBAgBMmDABERERrN2BhZUF+eHDh6iqqoKVlRVPj+/SpQu6dOlCm6rXIyAgAB4eHtzOeT/i6emJc+fO4fPnz0JOxj4/e7AA/N+5T+q/3rx5g7i4OO6t0j/SsmVLDB48GMeOHRNyMuFgZUGufcP/zP3qbP8qIyy8rAb4Nw0NDdja2uLEiRNCTMZOjx8/RklJCWxsbHh+zsiRI/HgwQNkZmYKLxhLBQYGYvz48T/VGIzNH3CsK8hfvnzBqVOnMHHixJ96nru7O+7cuYOsrCwhJWOny5cvQ1dXF8bGxj/1PPoBVzd/f/8frgb4t9o7Jf+5GzX19WAhICDgp5v329raIjc3F/Hx8UJKJjysK8inTp2CpaUldHR0fup5KioqcHd3x5EjR4SUjJ38/f0btFuFg4MDMjMzkZiYKIRU7FRWVoYTJ078cDVAXaZOnYqAgABWL9kStOvXr0NdXZ3b6J9XMjIymDJlCisPGFhXkH/2/Nw/0SVb3/r48SOioqIwduzYn36urKwsXbL1L2FhYejdu3eDevGamZlBWVmZuxs1xd/f+uTJkxESEsK6TSlYVZCfP3+OtLQ0bt/jn9WvXz/Iysrizp07Ak7GTkFBQXB2dkbTpk0b9PwpU6YgKCiILtn6/753a++P1C7Zoh9wXxUUFODSpUsYP358g57fvn17mJqacvt4swWrCnJAQAA8PT0hJyfXoOc3hnWKglK7GqChBQQAOnXqBENDQ1y4cEGAydjp5cuXSEhI4GsjhAkTJuDChQsoLCwUYDJ2Cg4OhpOTE5o3b97gMdi4PJM1BflHjYR45enpibNnz6KoqEhAydjp0aNHKC8vr7eREK/Y+KYXhsDAQEyYMKHORkK8UldXh4ODA204hIZf2/gnV1dXxMbG4u3btwJKJXysKciXLl1C27ZtYWhoyNc4mpqaGDRokMQv2WrI0sG6uLu7S3zDodqlg/x826hFGw4BT548QWFhIQYOHMjXOIqKihg7diyreq+wpiD/zN06PyLpR3WlpaU4efLkdxsJ8UpJSQmjR49m1Zte0K5evQptbW2BNJm3t7dHTk4Odz8+SdSQpYP1YdvqFVYU5JycHNy4ceOHjYR4NWTIELx9+xbJyckCGY9tTp06BXNzc+4mpvyq/YBjy5te0Pg9F/9Pkt5wqKysjKdGQrzq0aMHmjZtihs3bghkPGFjRUGubSSkqqoqkPFkZWUxadIkiX3TC+L83D/17t0bioqKEtlwKC8vD1evXm3Q0sH6TJ48GcHBwaioqBDYmGxx5swZmJmZoU2bNgIZj20X8sW+IAtiNUBdvLy8cPToUYlbsvXixQskJyfz1EiIV5LccCgoKAjDhw/nqZEQrzp06IDu3btLZMMhQZ6arDVhwgRcvHiRFQ2HxL4gP3jwADU1NbC0tBTouJ07d4a+vj7Cw8MFOq64+9lGQrzy8PDA+fPnJarhkLAOFgDJvM7x+vVrxMfH89xIiFctWrTAkCFDEBISItBxhUHsC7KgVgPURdLe9NXV1T/dSIhXGhoasLOzk6glWz/ag5Afbm5uePjwITIyMgQ+triqbSTEz9LB+rDltIVYF+SGNhLilbu7O2JiYiRmydaVK1fQunXrn24kxCu2vOkF5Ud7EPKjsewRx6uGNhLila2tLfLz8xEXFyeU8QVFrAvyyZMn0b9/f7Rq1Uoo4ysrK8Pd3V1i3vT89AbghSTtEcfLHoT8agx7xPHq2rVrP9yDkB/S0tKYMmWK2H8jFuuCLIwT/P9Wu06xsTcc+vjxI65duybQ1QD/JklLtnjdg5AfPXv2hKqqqkQ0HBLWufh/mjx5Mo4dOybWDYfEtiA/f/4cz58/h5OTk1Dn6dOnD+Tl5XH79m2hzsO0o0ePwsXFpcGNhHglKXvECXrpYF3YtmSroX5mD0J+tGvXDj169MDZs2eFOg8/xLYg+/v789VIiFeS8KYX5mqAf5OEPeJq9yDkp5EQryZMmIDw8PBG3XAoODiY5z0I+SXuF/LFsiBXV1fjyJEjIikgwNclW+fOnWu0DYcePnyIyspKnvcg5Je4v+n5FRAQgAkTJgh86WBd2L5H3I80ZA9Cfri6uuLJkydi23BILAtyZGQk2rdvD319fZHMp6mpCVtb20a7ZEuYSwfrMnLkSNy/f79R7hEnyEZCvGrM3+CePHmC4uLin9qDkB8KCgoYO3as2G6XJZYFWRTn5/6tsTYH//LlC06ePCm0pYN1UVJSarRLtq5cudKgPQj5weY94n5EkI2EeCXODYfEriDn5OTg5s2bGD16tEjndXBwQEZGBpKSkkQ6r7DV7kEoqEZCvGqsDYdE+fW6Vu0ecY3tgEHQjYR41aNHDzRv3hzXr18X6by8ELuCfPToUbi6ugqskRCvaveIa2xfDZn4tgEAvXr1anR7xPGzByG/aveIa0wNh/jZg5Bf4noaSKwKsihXA9Slse0Rl56ejtTU1AbvQciPxthwiN89CPnRvn17mJiYsG6PuO8RxX0G9Rk/fjwiIyNRUFDAyPz1EauCfO/ePRBCYGFhwcj8jW2PuNqlg6JYDVAXDw+PRrNHHNMHC0Dj2k3k1atXePbsmUiWDtalRYsWGDp0qNg1HBKrglz79VpUqwHq0liO6gS1ByE/1NXVYW9v3yhWrwhqD0J+1O4R9+7dO8YyCIog9iDklzheyBebglxSUoLTp0+LdDVAXRrLHnGC2oOQX+J6ru5niXrpYF1q94gT1yVbvKpdOsjkwQIADBo0CAUFBWLVcEhsCvLJkydhbW0NbW1tRnMoKSlh1KhRrF+yxdTFvH+zt7fHhw8fWL1HnCD3IOSXOC/Z4lVUVBS0tLTQvXt3RnPUNhwSpwMGsSnITCwnqg/bu2zl5OTg+vXrAtuDkB+NoeGQoPcg5Afb9oirizj9rU+ePBnHjx8Xm4ZDYlGQ09LS8OLFCzg6OjIdBcDXPeIUFBRY23Codg9CNTU1pqMA+L+GQ2xdsiUu3zaA/+u9wtYPuLy8PFy5ckXojYR41bZtW/Ts2RNnzpxhOgoAMSnI/v7+mDhxotAbCfGKzQ2HxGE1wL916NAB3bp1Y2XDoRcvXiAlJUWgexDya/z48YiIiGDFHnH/JspGQrwSpwv5jBfkqqoqHDlyRGyOQGqxdY+4+/fvo7q6Gv3792c6yjfY+gHn7+8vlD0I+dGyZUvW7BH3T+J4sAAAI0aMQFxcHN68ecN0FOYLcmRkJDp06CCyRkK8YuseceKwdLAubm5uePToEav2iBOHpYP1EaejOl49fvwYX758wYABA5iO8g0FBQWMGzdOLFavMF6QxfETsxbbjupKSkpw6tQpsVgN8G+1e8SJw5ueV5cvX4aenh6MjIyYjvIfdnZ2yMvLw9OnT5mOwrNDhw6JvJEQr2pXr9TU1DCag9FX5sOHD7h165bIGwnxim17xJ06dUqoexDyy8vLi1VLtsTpYt6/sWWPuFqlpaU4ceKEyBsJ8crU1BQtW7ZkvOEQowX5yJEjcHNzg4qKCpMx6sW2JVvi/G0DAMzMzKCqqoqbN28yHeWHcnNzhb4HIb9qGw6Jy5Kt7wkLC0Pfvn2Fugchv8ThGzFjBZkQwmhzEV6xZY+458+fIz09Xeh7EPKDTQ2HgoKC4OLiIjZLB+tSu0ccGxoOseFvffz48bh06RKjDYcYK8gxMTEAAHNzc6Yi8IQte8SJag9Cfnl4eIj9HnHiuhqgLmxoOPTy5UskJiaK1dLBujRv3hyOjo4IDg5mLANjBbn2E1PcVgPURdyP6mr3IBTX853/1LJlSzg4OIj1HnGi3oOQH+K+RxwgHo2EePXPv/Xq6moQQkQ6v0gLcl5eHubOnYvi4mKEhYXB09NTYGO7ubkhOzsbwNfVBvb29gIb+597xK1cuRLPnz8X2Nj8CA8Px/HjxxEZGYl27drBwMBAIOO+efMG48eP574Z79y5g8WLFwtkbOD/ztVVVVUx3kzqn/744w+8efNG4I2EAgICsH//fu7/r1mzBhcvXhTI2P/cI+7JkyfYsmWLQMblV2VlJSZPniyURkJTp05FcnIygK9F09HRUWAbFA8aNAiFhYV48uQJXFxcEBsbK5BxeSXSglxRUYETJ04gNDQUAwYMQHV1tcAWY7dp0wabNm0CABw8eBAtW7YUyLhZWVnIzMzk7hF39uxZlJSUCGRsfr1//x7Xr1/nfr2ubRHJr9atW+Px48eIi4sDIQS///47TExMBJAYiIuLQ58+fZCbm4vr168jKipKIOMKQkJCAh48eIBTp07Bw8MDd+7cEci4vXv3xvLly1FVVYWioiL4+fmhR48efI9LCMGdO3e4q1fu37+PtLQ0ASTmn4yMDIKDgxEZGYlWrVpBS0uLW0T51bVrV6xevRrA14OSkpISgeww9ODBAwQHB3NXr7x8+RLKysp8j/tTiAhVVlYSWVlZYm5uTvz8/EirVq3IuXPnBDJ2VlYWad68OVFSUiIaGhokMTFRIOPevXuXaGhokF27dpEOHToQLS0tkpWVJZCx+XX27Flib29PmjZtSv7++2+io6NDPnz4IJCxjx49Srp37066d+9OunTpQqqrqwUy7qpVq4iJiQmZN28eGTNmDOnZs6dAxhUEb29vMmnSJDJkyBDi5uZGnJ2dBTb26NGjiZWVFbGysiJz584VyJgcDof06tWLzJgxg5iYmBAPDw+yfPlygYwtCNra2sTJyYmsXr2adOjQgezfv18g4xYXFxNNTU3Spk0b0q5dO3Lt2jWBjPv27VvSsWNHMm/ePNK8eXPStGlTUlBQIJCxeSXSgkwIIS1atCBNmzYlGhoa5NSpUwIde86cOURKSoq4ubkJdNwbN24QDQ0NoqenR2RkZARWnPj14MEDoqurS0xNTUnHjh3Jy5cvBTZ2VVUV0dPTI02aNCFBQUECG5fD4ZBVq1aRtm3bEmVlZTJ06FCBjc2v5cuXE11dXWJsbExGjx5NysvLBTZ2QkICUVRUJIqKigL9QP/8+TMZOHAg6dGjB2nXrh3ZvXu3wMbmV7du3YiioiLR1tYm+/btE+jYGzduJLKysqRXr16Ew+EIbNwPHz4QU1NToqurS+Tk5AQ6Ni9EXpA1NDSIvLw8uX79usDHzsrKItLS0iQ2NlbgY8fFxRFVVVWioKAg8LEb6s2bNwQA6dy5s8COjP9p7dq1REVFRSgfQHv27CFSUlJk4JthpIMAABVJSURBVMCBAh+7oVauXEkAkP/9739C+TcbGxsTe3t7gY9bVlZGhg0bRgCQI0eOCHz8hjIwMCAyMjLk9OnTAh+7uLiYyMnJkRMnTgh87MLCQmJoaEgUFRUFPvaPiLwgT506lURGRgpt/C9fvght7MePH5MxY8YIbfyfVVlZSQYNGkQ+ffoklPE5HA4pLS0VytiEELJp0yayZcsWoY3/s65du0a8vLyEdlRUXl4utG9X1dXVxNXVVWCn6gTBx8eHHD58WGjjC/Nvvbi4mBw/flxo49dHihARr+ugKIqi6iQrjEEJIXj9+jXS09O5P1NXV0f37t2FcuMCh8NBamrqN53EtLS00K1bN8jIyAh8vurqaiQmJiInJ4f7Mz09Pejr6wulcUpVVRWePXuGvLw87s86d+6M9u3bC2Udd0VFBeLi4ritR6WkpGBgYAA9PT2BzwV8XaYYFxeH0tJSAF/7NBgbGwutJ0dhYSGePn3KbZgvJyeH7t27Q11dXSjz5eXlIT4+HtXV1QCAJk2acHf+EIacnBwkJCRwG+UoKyujZ8+eUFJSEvhchBBkZGQgNTWVu0xSTU0NPXr0gIKCgkDn+fjxI169elVnS9wmTZqgXbt2aN26NWRl+S9rHA4H2dnZeP36Nb58+fKf3ysqKqJdu3bQ1dUVaI0RaEEuLCzE/PnzER4eDmlpaXTo0IFboHJzc5GVlYU+ffrAz89PIMuo3r59izlz5iA6Ohpqampo3bo1pKSkQAhBTk4OcnJyYG5ujq1btwpks8/k5GTMnTsX9+7dg5aWFrS0tLjzZWZmoqioCAMGDMC2bdvQtm1bvueLj4/HvHnz8PDhQ+jo6EBTUxPA1zfLq1evwOFwMGzYMPz9998Cafh969Yt+Pj44OnTp2jfvj2aN28O4OumlM+fP4eioiLc3d2xfv16KCoq8j3f8ePHsWHDBjx//hydO3fmLl2qqqpCWloamjZtijFjxmDNmjV831RACMGePXuwc+dOvH37Fl27duUWqIqKCqSlpUFTUxNTpkyBj48P339kNTU12LBhAwICApCbm4uuXbty/w2lpaVIS0tD27ZtMXPmTMyYMYPvD9bq6mqsXr0aQUFByM/Ph76+Pvfgp6ioCC9evICBgQFWrlwpkDvmiouLsWjRIpw7dw7V1dXo1KkT9zUrLCzEq1evYGxsjDVr1mDw4MENmuPdu3dYsWIFHj58iDdv3qBJkybQ1dWFmpraf16viooKZGVlIS8vD7q6uujcuTN+//132NjY8DxffHw8Vq1ahYSEBGRmZkJFRQWtW7euc+lbWVkZ3r9/j8LCQujq6qJr165YunQpLC0tG/RvrSXQUxbOzs6QkpKCt7c3tLW1//OiFRcXIyoqCgcOHEBKSgpatGjR4LmqqqrQvXt3DBgwAG5ubnWuO/78+TMuX76Mo0ePIjU1la++BEVFRdDX14enpycGDx5c59FNfn4+wsLCEB0djWfPnvH1bSA/Px+Ghobw9vaGnZ3df9ZZEkLw4cMHHDhwAIQQvm/tfv36NczMzLB48WJYWlr+52iKEIJ3795h+/bt6Nq1K/bt28fXfDdu3MC4ceOwYsUKmJqa/qcBPIfDwdu3b+Hn58f9EOdHcHAwli1bhmXLlsHQ0PA/R1E1NTVIT0/Hxo0bMXHiRCxcuJCv+Xx9fXHkyBH4+Ph8U6xqVVdXIzk5GWvWrMGaNWswYcIEvuZbvXo1zp07h8WLF39zIFSrvLwcsbGxWLNmDS5fvgwzMzO+5hszZgyKiorwyy+/QFdX9z9/62VlZXjw4AH++usvREdHw9jY+KfGr6iogLGxMaysrGBnZwddXV2empBVVlYiOzsbiYmJ2L59O65evcrTmu/s7Gx069YNXl5e6NWrF3R0dHg66KioqEB2djaePXuGHTt24O7du3zdoCWwgvz582fo6OggKirqh7srLF68GF5eXny9Ce/du4cpU6YgKCjoh0cXs2fPxsKFC+Hq6trg+cLCwvD3339j27Zt330cIQSenp4ICAhAv379GjxfUFAQAgMDsXHjxu8+rrKyEnZ2dsjKyuLrK/DmzZsRGxv7wzvyPn36hBEjRqC4uJivo7pJkyZBW1v7h61Xs7OzMXny5G9O1zSEnZ0dhgwZgoEDB373cc+ePcOWLVuQmJjI13xGRkZYtGgRunXr9t3HXb9+HVeuXMHVq1f5mq9Dhw7466+/0Llz5+8+7uDBg2jSpAm2bt3a4LnKy8vRvHlzREVF/fC0xK5du9CqVSusW7fup+aIjo7GzJkz+eqffeDAAaioqMDX1/eHj92zZw8iIyOxYsWKBs+3Y8cOtGnThnvTSkMI7IRnSkoKOnTowNNWNx07duT7rp2UlBR07tyZp6LQqVMnvudLTk5Gp06dfvg4KSkpdOrUCSkpKXzP17Fjxx8+Tl5eHu3bt0dqaipf8yUl/b/2rjWmzfJ9X9giXXCZlnKYUGzpyqlSThtOVmSMsYHCQCSgSOY3SZhbMAaMW5ZtLhqPX9REl8WoMxTHJguKbBB24FDOuLVSNlihhQkjlVI2EMra0t+Hpu9/Ly2j7fti/P9+vRI+9OXlvbjfPvdzuJ/7vh7lms4MWAVYfH198eeff1LiGxwcRHh4+Jr3BQUFwWg0Uu6Qb9265ZR9QqEQIyMjlITKzWYzRkdHnWov4eHhlL87g8GAu3fvgs/nO8VH1ReGh4cREhLiVIw4PDwcSqXSZQ6bf1OBUCh0mttZf6OLbzXQ2iHzeDzStY6ODuTn5yMvL4800vH5fMr/+ODgIClOe+LECWRkZDiccfF4PMp8SqWSZN/U1BRKS0tRUFCAwsJCklhOaGgo5UavVCpJDra0tIT9+/fjtddeQ2FhISlkwOPxaBkAHDm02WxGcXExysvLiWthYWGU+CwWC4aGhuzaS05ODoqKilBcXEzonHh5eVHmm5+fh06ns9sknJubQ2VlJV555RUUFBRAoVBgw4YNYLPZlEr61Wo1/Pz8SEtemz6I7Sc1NRVSqRSbN2/G9PQ0pXL84eFhBAcHk8IwVVVVKCwsRGFhIQ4fPkxsYNIxeN+8edOurVRXVxN8D5/1x+fz3frunPHve/fuoaysDC+//DLKysrs9Cx4PJ7Ttq70N0d8zc3NKCwsxLZt2xz6Nx1+SFuHPDQ0RBKfNpvN+Pjjj/HFF1/g3LlzaGxsxOjoKADr0dtUBXqGhoYQGhpKfM7JycGXX37p8F4ej0eZb3h4mNRAmEwm3n77bZw/fx7fffcdzp07R7KPaqO/ffs2yb7HH38c33zzDaqrqyGVStHR0UGcZMLlcilrGKhUKhKfDdXV1XbOx+VyKb1PrVYLJpPpcCPy1KlTkEql+PHHH4lroaGhlPhstq2M43722WdITk7Gzz//TLKTz+dT4hseHrYbbHg8HqRSKWEbi8VCWloaGAwGQkNDoVKp3Oa7ffs2qW1qtVqcPXsWZ86cQU1NDZaXl9HU1AQAePrpp6HVailpnqz0dZVKhQsXLuDMmTOQSqVob2/H+Pg4AGtb0Wg0Lp8Sc+vWLZJNjvz7+++/R1JSEi5cuICkpCS78EZISAgmJyed0jJXqVRr8gkEAnzyySerxqR5PB5GRkYoKcTR1iEbDAbSEkapVILL5SIkJATe3t7Ys2cPWlpaAFgVqqgKvi8tLZH4EhISVt20Y7FYlEV3VvJxOBziYFZfX1/weDxotVqCzzYjcRcPHjwg8Xl5eREbbSaTCSaTiQjXsFgsLC4uUuJbaR9gTZ+SyWTIy8sjXaf6PpeWllzK0vDx8aH0Ph3ZZku1y83NBWBNfbNtnNJh36OW8729vQgODiZm7FTby9LSkl0WitlsxtLSEkwmEwwGA/z9/QFYRX+8vb0p+d9KPo1Gg5iYGLBYLDCZTCQkJODq1asArN+d2Wx2uUN2xr9bWlqQnZ0NAMjOzrY7iYbJZOKxxx4j0g2p8vH5fLuB9mHQ0a+tm9qbVqtFYGAg8TkgIIDosP7bMDk5iaGhIZd3kl2FLXyQkZGB5557bt35Pv/8cxw6dOgf06z28vLCgQMHUFJSgtra2nXlmpiYwJNPPokTJ06guLgYJ0+epDyoOYvGxka3U8GcQUBAAEpKSpCdnY3MzEw88cQTlDaY14JAIMD169cxOzsLg8EAmUxGytFfL8zMzBC54xwOB3q9ft051xv/qPzm/wcxelexsLCAyspKvPPOO+t+NiCDwYBUKkVDQwOUSiWlZe5aaGtrA5vNpk1j2Rl8++23qKqqIsJcv//++7pxmc1mDA0NoaCgAFKpFBs2bPhHTsQ2Go1obW3F7t27143j/v37aGlpwS+//IJLly5hcXGRNv1lR+Dz+di/fz8OHDiAgwcPQigUrktB1v8C1q1DDggIII2SWq2WWDbRlfrs7HPWi89kMqGyshKZmZnYtWvXuvPZsHHjRiQmJqKzs5MWHkd8crkcra2tyMnJwZEjR9Db24ujR48+8v9yl8sGW/tgs9nYuXMn5Y3YR/EFBAQgICCAWGWkp6cTcf/1sg8AZDIZIiMjSXnzdKsX2AqJnnrqKTCZTKSlpUGhUKwbHwDk5eWhqqoKp0+fxqZNm0j7Ee7yrfV3bDabyL6Znp4mCplceYa7967H3wM0dsj+/v6kwwGjo6Nx584dTExMwGg0oqmpCS+88AIAa9EDVQF5Pz8/pw8j1Ol0hLO7Cw6HQ+KzWCx4//33wefzUVJSYsdHtQyXzWaT+PR6Pebm5gBY4/U9PT1EPEun0xFVfO5i5ft866230NDQgF9//RUffPABtm3bhpMnTwKwLhWpvE8/Pz/o9XpSatni4iJRomorKrClIVF9n35+ftDpdKRrHA4HgYGBRDZFT08PwsLCCD4q9q1sKw/DUbiCDvse5gsKCsLAwAAMBgMsFgt6e3uJtrKwsACz2UxJeN3f398uPGDjn5qawpUrVwgb9Xo9Nm7c6PKM+VHv0IbU1FTU19cDsArVp6amkn4/Pz8PBoPhVHoem822ayOuQqfTgc1mU4oE0FY6HR0dTUpuZzKZqKiowMGDB2E2m7Fv3z7CwdRqNeVSZpFIRNLKOHz4MPr7+zE7O4sXX3wRb775JrEZRQdfdHQ01Go1cSirXC5HQ0MDtmzZguLiYgBAWVkZJBIJNBoNRCIRLXwJCQkArDOAY8eOYXl5GcvLy8jIyCDOfBsfH0d+fj4lvsjISKjVaqf0KjQaDaVQhq+vLzgcDu7evUvs1ut0OlRUVACwhhP27t2L5ORkWvj4fD5mZmawsLBAqkCsqKjA0aNHYTQaERwcjGPHjsFisWB0dJQSX1RUFEZHR2GxWEjOaRtIjxw5QlxbWFjAzMzMIzeLnOWz4dlnn0V6ejpef/11MBgMREREEO1Do9E4rBx0le/s2bOka5WVlbh37x6YTCbeffddYkNMrVYjIiLC5U5KJBJBrVYTnx359xtvvIH33nsPdXV1CAoKwkcffUR6hlqthlAodEpfJioqChqNhtiod8S3adMmfPrpp9Dr9SgvL0d4eDi++uor4hkajQYREREu2bkStHXIUVFRpBcIABKJBBKJxO7esbExoqNxF9HR0aTjfz788MNV7x0bG0N6ejolPpFIROwcA0BcXNyq522Nj49Tqgq08V2/fp34LBQKSfmdD0OtVlOO9doGANsq5mFs3boVW7duBWAtaVar1UTDdRe2AcDWIYeEhDg8+NRoNGJiYsKpIpLVwGAwEBYWhrGxMdJ7ioiIIKXXAcBff/0FHx8fSis4DocDb29vu5kvi8XC5cuXSfdqNBoIBAJKHSSXy8Xff/+Nubk5IlOktLQUpaWldvfS0VZWDgCAtQLQEdydDEVHR6Orq4v4vJp/f/3116s+wxVbVw4Aq/E9qtKTjokfbSELgUCA2dnZNRPqTSYT+vv7KZ8pFhcXh4GBgTUT6o1GIy188fHx6Ovrg9FofOR98/PzGBgYoCyeFB8fj/7+/jVTdjQaDWZnZ4nltrtITExEb2/vmnEwuVyOwMBASroggDWt6GGHWw19fX0QCASUxYWc5evs7KRF+Co2NhYdHR1r3tfd3U15cuLl5QWxWOzUnkJ3dzdlHQsulwuTyeTUpnJPT49b9sXFxeHGjRsOldacgcViccnWhIQE9Pb2ul2h6SrfaqBVXOjUqVM4fvw4srKyIBaL7dTe5HI5rl27hmeeeQZ1dXWUsy7Ky8tRX1+P3bt3QywWg8vlEuprU1NTUCgUuHr1KkQiEX766SdKfBaLBUVFRRgcHERaWhrEYjEhoGSTIFQoFGhubkZ2djYlrQAb3759+3Dnzh2kpqYiNjbWTu1NoVDg4sWLOH78uMPZkCswGo1ISUmBl5cXduzYAbFYTIg/2dTe/vjjD1y8eBGnT5+2y012FdPT00hMTERkZCQSExMRExNDUnsbHh6GQqHApUuXUFVVhczMTEp8IyMjSE5Oxvbt2xEXFweRSERSe7t58ybkcjmuXLmCxsZGYkXgLvr6+rB3717s2rULsbGxiIqKIqm9KZVK3LhxA11dXejo6KBcttve3o7c3FxkZGRALBYjMjKSkDG4f/8+cYDrxMQE+vr6KB8K+sMPP6CyshJZWVmIiYmxU3tTKBTo6urCzMwMenp63JL+PHToEGpra/H8889j8+bNCA4ORkhIiEO1N4PBgMnJSUxMTGBychKDg4N48OABOjs7ncp+MpvNyMrKwvj4OOLj4wm+1USNFhYWSHwDAwNgMpmQyWSUZEdpF6jv7u7Gb7/9hvb2dlKM18/PDzt27MDOnTuRn59PS1qMxWJBU1MTLl++jPb2djs9ZBtfbm4uLTrFy8vLqKurw7Vr1+xyLblcLiQSCdLT07Fnzx5aUvzMZjNqa2sJvoc3HYRCISQSCbKzs5GUlESZC7B2TDU1NWhtbUVHRwdRimrTQ05JSUFubu6agjnOQq/Xo6amBjKZDN3d3YQeMoPBgEgkgkQiQV5eHm2pd1NTUzh//jza2trQ399P0kOOjY2FRCJBfn6+U5oQzkCtVqO2thbt7e2Qy+XE6srHxweJiYlISUlBQUEBgoKCaOFTqVSora1FW1ubnR7y9u3bIZFI8Oqrr9KWntnf34/6+nrIZDI7PWQbX1FRESUd5u7ubvT19UGlUmFkZAQjIyN2JdLA/+khCwQCbNmyBeHh4XjppZdc0ka2WCxobm7GwMAAwbWaHjKLxQKfz4dAIIBAIEBkZCSysrIo92ueE0M88MADD/4l+EcLQzzwwAMPPFgdng7ZAw888OBfAk+H7IEHHnjwL8F/AEvybU5yY8XgAAAAAElFTkSuQmCC">
</div>

</div>

</div>
</div>

</div>
<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>Now that we finished setting up our model we can start calculating the formulae of interest. Let's first review the formulae for the share and elasticities in the Nested Logit Model. The share of one end node in the Nested Logit Model is
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-1-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mi&gt;j&lt;/mi&gt;&lt;/msub&gt;&lt;mo&gt;=&lt;/mo&gt;&lt;mfrac&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;4&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;5&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;6&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;7&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;8&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;10&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;11&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;9&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;4&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;5&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mrow&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;4&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;5&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;6&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;7&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;8&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;10&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;11&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;9&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;4&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;5&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;6&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;7&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;8&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;10&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;11&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;9&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;4&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;5&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;msup&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;v&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;/msup&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;.&lt;/mo&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-1" style="width: 69.586em; display: inline-block;"><span style="display: inline-block; position: relative; width: 57.979em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(-8.092em, 1057.92em, 17.324em, -999.997em); top: -2.199em; left: 0em;"><span class="mrow" id="MathJax-Span-2"><span class="msubsup" id="MathJax-Span-3"><span style="display: inline-block; position: relative; width: 0.836em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-4" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="mi" id="MathJax-Span-5" style="font-size: 70.7%; font-family: MathJax_Math-italic;">j</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-6" style="font-family: MathJax_Main; padding-left: 0.301em;">=</span><span class="mfrac" id="MathJax-Span-7" style="padding-left: 0.301em;"><span style="display: inline-block; position: relative; width: 55.241em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(7.622em, 1053.28em, 17.622em, -999.997em); top: -17.914em; left: 50%; margin-left: -26.664em;"><span class="mrow" id="MathJax-Span-8"><span style="display: inline-block; position: relative; width: 53.277em; height: 0px;"><span style="position: absolute; clip: rect(2.979em, 1053.28em, 8.336em, -999.997em); top: -6.009em; left: 0em;"><span class="msubsup" id="MathJax-Span-9"><span style="display: inline-block; position: relative; width: 53.277em; height: 0px;"><span style="position: absolute; clip: rect(2.979em, 1052.03em, 8.098em, -999.997em); top: -5.771em; left: 0em;"><span class="mrow" id="MathJax-Span-10"><span class="mo" id="MathJax-Span-11" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-12"><span class="msubsup" id="MathJax-Span-13"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-14"><span class="mo" id="MathJax-Span-15" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-16"><span class="msubsup" id="MathJax-Span-17"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-18"><span class="mo" id="MathJax-Span-19" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-20"><span class="msubsup" id="MathJax-Span-21"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-22" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-23"><span class="mrow" id="MathJax-Span-24"><span class="mfrac" id="MathJax-Span-25"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-26"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-27" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-28"><span class="mrow" id="MathJax-Span-29"><span class="mn" id="MathJax-Span-30" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-31"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-32" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-33"><span class="mrow" id="MathJax-Span-34"><span class="mn" id="MathJax-Span-35" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-36" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-37" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-38" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-39" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-40"><span class="mrow" id="MathJax-Span-41"><span class="mfrac" id="MathJax-Span-42"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-43"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-44" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-45"><span class="mrow" id="MathJax-Span-46"><span class="mn" id="MathJax-Span-47" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-48"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-49" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-50"><span class="mrow" id="MathJax-Span-51"><span class="mn" id="MathJax-Span-52" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-53" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-54" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-55" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-56" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-57"><span class="mrow" id="MathJax-Span-58"><span class="mfrac" id="MathJax-Span-59"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-60"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-61" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-62"><span class="mrow" id="MathJax-Span-63"><span class="mn" id="MathJax-Span-64" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-65"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-66" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-67"><span class="mrow" id="MathJax-Span-68"><span class="mn" id="MathJax-Span-69" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-70" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-71" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-72"><span class="mrow" id="MathJax-Span-73"><span class="mfrac" id="MathJax-Span-74"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-75"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-76" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-77"><span class="mrow" id="MathJax-Span-78"><span class="mn" id="MathJax-Span-79" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-80" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-81"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-82" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-83"><span class="mrow" id="MathJax-Span-84"><span class="mn" id="MathJax-Span-85" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-86" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-87" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-88" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-89"><span class="mo" id="MathJax-Span-90" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-91"><span class="msubsup" id="MathJax-Span-92"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-93" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-94"><span class="mrow" id="MathJax-Span-95"><span class="mfrac" id="MathJax-Span-96"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-97"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-98" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-99"><span class="mrow" id="MathJax-Span-100"><span class="mn" id="MathJax-Span-101" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-102"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-103" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-104"><span class="mrow" id="MathJax-Span-105"><span class="mn" id="MathJax-Span-106" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-107" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-108" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-109" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-110" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-111"><span class="mrow" id="MathJax-Span-112"><span class="mfrac" id="MathJax-Span-113"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-114"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-115" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-116"><span class="mrow" id="MathJax-Span-117"><span class="mn" id="MathJax-Span-118" style="font-size: 50%; font-family: MathJax_Main;">4</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-119"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-120" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-121"><span class="mrow" id="MathJax-Span-122"><span class="mn" id="MathJax-Span-123" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-124" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-125" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-126" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-127" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-128"><span class="mrow" id="MathJax-Span-129"><span class="mfrac" id="MathJax-Span-130"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-131"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-132" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-133"><span class="mrow" id="MathJax-Span-134"><span class="mn" id="MathJax-Span-135" style="font-size: 50%; font-family: MathJax_Main;">5</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-136"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-137" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-138"><span class="mrow" id="MathJax-Span-139"><span class="mn" id="MathJax-Span-140" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-141" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-142" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-143"><span class="mrow" id="MathJax-Span-144"><span class="mfrac" id="MathJax-Span-145"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-146"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-147" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-148"><span class="mrow" id="MathJax-Span-149"><span class="mn" id="MathJax-Span-150" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-151" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-152"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-153" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-154"><span class="mrow" id="MathJax-Span-155"><span class="mn" id="MathJax-Span-156" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-157" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-158" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-159"><span class="mrow" id="MathJax-Span-160"><span class="mfrac" id="MathJax-Span-161"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-162"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-163" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-164"><span class="mrow" id="MathJax-Span-165"><span class="mn" id="MathJax-Span-166" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-167" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-168"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-169" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-170"><span class="mrow" id="MathJax-Span-171"><span class="mn" id="MathJax-Span-172" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-173" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-174" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-175" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-176"><span class="mo" id="MathJax-Span-177" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-178"><span class="msubsup" id="MathJax-Span-179"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-180"><span class="mo" id="MathJax-Span-181" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-182"><span class="msubsup" id="MathJax-Span-183"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-184" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-185"><span class="mrow" id="MathJax-Span-186"><span class="mfrac" id="MathJax-Span-187"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-188"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-189" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-190"><span class="mrow" id="MathJax-Span-191"><span class="mn" id="MathJax-Span-192" style="font-size: 50%; font-family: MathJax_Main;">6</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-193"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-194" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-195"><span class="mrow" id="MathJax-Span-196"><span class="mn" id="MathJax-Span-197" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-198" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-199" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-200" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-201" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-202"><span class="mrow" id="MathJax-Span-203"><span class="mfrac" id="MathJax-Span-204"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-205"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-206" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-207"><span class="mrow" id="MathJax-Span-208"><span class="mn" id="MathJax-Span-209" style="font-size: 50%; font-family: MathJax_Main;">7</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-210"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-211" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-212"><span class="mrow" id="MathJax-Span-213"><span class="mn" id="MathJax-Span-214" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-215" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-216" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-217" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-218" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-219"><span class="mrow" id="MathJax-Span-220"><span class="mfrac" id="MathJax-Span-221"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-222"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-223" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-224"><span class="mrow" id="MathJax-Span-225"><span class="mn" id="MathJax-Span-226" style="font-size: 50%; font-family: MathJax_Main;">8</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-227"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-228" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-229"><span class="mrow" id="MathJax-Span-230"><span class="mn" id="MathJax-Span-231" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-232" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-233" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-234"><span class="mrow" id="MathJax-Span-235"><span class="mfrac" id="MathJax-Span-236"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-237"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-238" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-239"><span class="mrow" id="MathJax-Span-240"><span class="mn" id="MathJax-Span-241" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-242" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-243"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-244" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-245"><span class="mrow" id="MathJax-Span-246"><span class="mn" id="MathJax-Span-247" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-248" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-249" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-250" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-251"><span class="mo" id="MathJax-Span-252" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-253"><span class="msubsup" id="MathJax-Span-254"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-255" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-256"><span class="mrow" id="MathJax-Span-257"><span class="mfrac" id="MathJax-Span-258"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-259"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-260" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-261"><span class="mrow" id="MathJax-Span-262"><span class="mn" id="MathJax-Span-263" style="font-size: 50%; font-family: MathJax_Main;">10</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-264"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-265" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-266"><span class="mrow" id="MathJax-Span-267"><span class="mn" id="MathJax-Span-268" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-269" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-270" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-271" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-272" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-273"><span class="mrow" id="MathJax-Span-274"><span class="mfrac" id="MathJax-Span-275"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-276"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-277" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-278"><span class="mrow" id="MathJax-Span-279"><span class="mn" id="MathJax-Span-280" style="font-size: 50%; font-family: MathJax_Main;">11</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-281"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-282" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-283"><span class="mrow" id="MathJax-Span-284"><span class="mn" id="MathJax-Span-285" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-286" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-287" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-288" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-289" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-290"><span class="mrow" id="MathJax-Span-291"><span class="mfrac" id="MathJax-Span-292"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-293"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-294" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-295"><span class="mrow" id="MathJax-Span-296"><span class="mn" id="MathJax-Span-297" style="font-size: 50%; font-family: MathJax_Main;">9</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-298"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-299" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-300"><span class="mrow" id="MathJax-Span-301"><span class="mn" id="MathJax-Span-302" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-303" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-304" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 9.051em;"><span class="texatom" id="MathJax-Span-305"><span class="mrow" id="MathJax-Span-306"><span class="mfrac" id="MathJax-Span-307"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-308"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-309" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-310"><span class="mrow" id="MathJax-Span-311"><span class="mn" id="MathJax-Span-312" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-313" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-314"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-315" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-316"><span class="mrow" id="MathJax-Span-317"><span class="mn" id="MathJax-Span-318" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-319" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-320" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-321"><span class="mrow" id="MathJax-Span-322"><span class="mfrac" id="MathJax-Span-323"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-324"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-325" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-326"><span class="mrow" id="MathJax-Span-327"><span class="mn" id="MathJax-Span-328" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-329" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-330"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-331" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-332"><span class="mrow" id="MathJax-Span-333"><span class="mn" id="MathJax-Span-334" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-335" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-336" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 5.777em;"></span></span><span style="position: absolute; top: -6.366em; left: 52.324em;"><span class="texatom" id="MathJax-Span-337"><span class="mrow" id="MathJax-Span-338"><span class="msubsup" id="MathJax-Span-339"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.336em, 1000.3em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-340" style="font-size: 70.7%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.866em; left: 0.36em;"><span class="texatom" id="MathJax-Span-341"><span class="mrow" id="MathJax-Span-342"><span class="mn" id="MathJax-Span-343" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-344" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-345"><span class="mrow" id="MathJax-Span-346"></span></span><span style="display: inline-block; width: 0px; height: 6.015em;"></span></span><span style="position: absolute; clip: rect(2.384em, 1048.1em, 6.908em, -999.997em); top: 0.063em; left: 0em;"><span class="msubsup" id="MathJax-Span-345-MathJax-Continue-1"><span class="mrow" id="MathJax-Span-346-MathJax-Continue-1"><span class="mo" id="MathJax-Span-347" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-348"><span class="msubsup" id="MathJax-Span-349"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-350"><span class="mo" id="MathJax-Span-351" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-352"><span class="msubsup" id="MathJax-Span-353"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-354" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-355"><span class="mrow" id="MathJax-Span-356"><span class="mfrac" id="MathJax-Span-357"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-358"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-359" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-360"><span class="mrow" id="MathJax-Span-361"><span class="mn" id="MathJax-Span-362" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-363"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-364" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-365"><span class="mrow" id="MathJax-Span-366"><span class="mn" id="MathJax-Span-367" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-368" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-369" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-370" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-371" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-372"><span class="mrow" id="MathJax-Span-373"><span class="mfrac" id="MathJax-Span-374"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-375"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-376" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-377"><span class="mrow" id="MathJax-Span-378"><span class="mn" id="MathJax-Span-379" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-380"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-381" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-382"><span class="mrow" id="MathJax-Span-383"><span class="mn" id="MathJax-Span-384" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-385" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-386" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-387" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-388" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-389"><span class="mrow" id="MathJax-Span-390"><span class="mfrac" id="MathJax-Span-391"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-392"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-393" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-394"><span class="mrow" id="MathJax-Span-395"><span class="mn" id="MathJax-Span-396" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-397"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-398" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-399"><span class="mrow" id="MathJax-Span-400"><span class="mn" id="MathJax-Span-401" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-402" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-403" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-404"><span class="mrow" id="MathJax-Span-405"><span class="mfrac" id="MathJax-Span-406"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-407"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-408" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-409"><span class="mrow" id="MathJax-Span-410"><span class="mn" id="MathJax-Span-411" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-412" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-413"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-414" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-415"><span class="mrow" id="MathJax-Span-416"><span class="mn" id="MathJax-Span-417" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-418" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-419" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-420" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-421"><span class="mo" id="MathJax-Span-422" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-423"><span class="msubsup" id="MathJax-Span-424"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-425" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-426"><span class="mrow" id="MathJax-Span-427"><span class="mfrac" id="MathJax-Span-428"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-429"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-430" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-431"><span class="mrow" id="MathJax-Span-432"><span class="mn" id="MathJax-Span-433" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-434"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-435" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-436"><span class="mrow" id="MathJax-Span-437"><span class="mn" id="MathJax-Span-438" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-439" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-440" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-441" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-442" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-443"><span class="mrow" id="MathJax-Span-444"><span class="mfrac" id="MathJax-Span-445"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-446"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-447" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-448"><span class="mrow" id="MathJax-Span-449"><span class="mn" id="MathJax-Span-450" style="font-size: 50%; font-family: MathJax_Main;">4</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-451"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-452" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-453"><span class="mrow" id="MathJax-Span-454"><span class="mn" id="MathJax-Span-455" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-456" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-457" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-458" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-459" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-460"><span class="mrow" id="MathJax-Span-461"><span class="mfrac" id="MathJax-Span-462"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-463"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-464" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-465"><span class="mrow" id="MathJax-Span-466"><span class="mn" id="MathJax-Span-467" style="font-size: 50%; font-family: MathJax_Main;">5</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-468"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-469" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-470"><span class="mrow" id="MathJax-Span-471"><span class="mn" id="MathJax-Span-472" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-473" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-474" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-475"><span class="mrow" id="MathJax-Span-476"><span class="mfrac" id="MathJax-Span-477"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-478"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-479" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-480"><span class="mrow" id="MathJax-Span-481"><span class="mn" id="MathJax-Span-482" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-483" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-484"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-485" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-486"><span class="mrow" id="MathJax-Span-487"><span class="mn" id="MathJax-Span-488" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-489" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-490" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; position: relative; width: 1.193em; height: 0px;"><span style="position: absolute; top: -5.771em; left: 0.003em;"><span class="texatom" id="MathJax-Span-491"><span class="mrow" id="MathJax-Span-492"><span class="mfrac" id="MathJax-Span-493"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-494"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-495" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-496"><span class="mrow" id="MathJax-Span-497"><span class="mn" id="MathJax-Span-498" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-499" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-500"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-501" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-502"><span class="mrow" id="MathJax-Span-503"><span class="mn" id="MathJax-Span-504" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-505" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-506"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-507"><span class="mo" id="MathJax-Span-508" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-509"><span class="msubsup" id="MathJax-Span-510"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-511" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-512"><span class="mrow" id="MathJax-Span-513"><span class="mfrac" id="MathJax-Span-514"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-515"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-516" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-517"><span class="mrow" id="MathJax-Span-518"><span class="mn" id="MathJax-Span-519" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-520"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-521" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-522"><span class="mrow" id="MathJax-Span-523"><span class="mn" id="MathJax-Span-524" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-525" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-526" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-527" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-528" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-529"><span class="mrow" id="MathJax-Span-530"><span class="mfrac" id="MathJax-Span-531"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-532"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-533" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-534"><span class="mrow" id="MathJax-Span-535"><span class="mn" id="MathJax-Span-536" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-537"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-538" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-539"><span class="mrow" id="MathJax-Span-540"><span class="mn" id="MathJax-Span-541" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-542" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-543" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-544" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-545" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-546"><span class="mrow" id="MathJax-Span-547"><span class="mfrac" id="MathJax-Span-548"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-549"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-550" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-551"><span class="mrow" id="MathJax-Span-552"><span class="mn" id="MathJax-Span-553" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-554"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-555" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-556"><span class="mrow" id="MathJax-Span-557"><span class="mn" id="MathJax-Span-558" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-559" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-560" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-561"><span class="mrow" id="MathJax-Span-562"><span class="mfrac" id="MathJax-Span-563"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-564"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-565" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-566"><span class="mrow" id="MathJax-Span-567"><span class="mn" id="MathJax-Span-568" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-569" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-570"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-571" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-572"><span class="mrow" id="MathJax-Span-573"><span class="mn" id="MathJax-Span-574" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-575" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-576"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-577" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-578"><span class="mrow" id="MathJax-Span-579"><span class="mfrac" id="MathJax-Span-580"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-581"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-582" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-583"><span class="mrow" id="MathJax-Span-584"><span class="mn" id="MathJax-Span-585" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-586"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-587" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-588"><span class="mrow" id="MathJax-Span-589"><span class="mn" id="MathJax-Span-590" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-591" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 5.182em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 10.658em;"></span></span><span style="position: absolute; clip: rect(12.92em, 1054.65em, 28.217em, -999.997em); top: -13.092em; left: 50%; margin-left: -27.557em;"><span class="mrow" id="MathJax-Span-592"><span style="display: inline-block; position: relative; width: 55.122em; height: 0px;"><span style="position: absolute; clip: rect(3.217em, 1054.17em, 8.813em, -999.997em); top: -6.247em; left: 0em;"><span class="mrow" id="MathJax-Span-593"><span class="mo" id="MathJax-Span-594" style="vertical-align: 2.86em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.598em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.604em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.068em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.473em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-595"><span class="msubsup" id="MathJax-Span-596"><span style="display: inline-block; position: relative; width: 53.277em; height: 0px;"><span style="position: absolute; clip: rect(2.979em, 1052.03em, 8.098em, -999.997em); top: -5.771em; left: 0em;"><span class="mrow" id="MathJax-Span-597"><span class="mo" id="MathJax-Span-598" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-599"><span class="msubsup" id="MathJax-Span-600"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-601"><span class="mo" id="MathJax-Span-602" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-603"><span class="msubsup" id="MathJax-Span-604"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-605"><span class="mo" id="MathJax-Span-606" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-607"><span class="msubsup" id="MathJax-Span-608"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-609" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-610"><span class="mrow" id="MathJax-Span-611"><span class="mfrac" id="MathJax-Span-612"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-613"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-614" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-615"><span class="mrow" id="MathJax-Span-616"><span class="mn" id="MathJax-Span-617" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-618"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-619" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-620"><span class="mrow" id="MathJax-Span-621"><span class="mn" id="MathJax-Span-622" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-623" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-624" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-625" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-626" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-627"><span class="mrow" id="MathJax-Span-628"><span class="mfrac" id="MathJax-Span-629"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-630"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-631" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-632"><span class="mrow" id="MathJax-Span-633"><span class="mn" id="MathJax-Span-634" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-635"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-636" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-637"><span class="mrow" id="MathJax-Span-638"><span class="mn" id="MathJax-Span-639" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-640" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-641" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-642" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-643" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-644"><span class="mrow" id="MathJax-Span-645"><span class="mfrac" id="MathJax-Span-646"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-647"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-648" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-649"><span class="mrow" id="MathJax-Span-650"><span class="mn" id="MathJax-Span-651" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-652"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-653" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-654"><span class="mrow" id="MathJax-Span-655"><span class="mn" id="MathJax-Span-656" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-657" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-658" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-659"><span class="mrow" id="MathJax-Span-660"><span class="mfrac" id="MathJax-Span-661"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-662"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-663" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-664"><span class="mrow" id="MathJax-Span-665"><span class="mn" id="MathJax-Span-666" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-667" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-668"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-669" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-670"><span class="mrow" id="MathJax-Span-671"><span class="mn" id="MathJax-Span-672" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-673" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-674" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-675" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-676"><span class="mo" id="MathJax-Span-677" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-678"><span class="msubsup" id="MathJax-Span-679"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-680" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-681"><span class="mrow" id="MathJax-Span-682"><span class="mfrac" id="MathJax-Span-683"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-684"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-685" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-686"><span class="mrow" id="MathJax-Span-687"><span class="mn" id="MathJax-Span-688" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-689"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-690" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-691"><span class="mrow" id="MathJax-Span-692"><span class="mn" id="MathJax-Span-693" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-694" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-695" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-696" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-697" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-698"><span class="mrow" id="MathJax-Span-699"><span class="mfrac" id="MathJax-Span-700"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-701"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-702" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-703"><span class="mrow" id="MathJax-Span-704"><span class="mn" id="MathJax-Span-705" style="font-size: 50%; font-family: MathJax_Main;">4</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-706"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-707" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-708"><span class="mrow" id="MathJax-Span-709"><span class="mn" id="MathJax-Span-710" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-711" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-712" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-713" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-714" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-715"><span class="mrow" id="MathJax-Span-716"><span class="mfrac" id="MathJax-Span-717"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-718"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-719" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-720"><span class="mrow" id="MathJax-Span-721"><span class="mn" id="MathJax-Span-722" style="font-size: 50%; font-family: MathJax_Main;">5</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-723"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-724" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-725"><span class="mrow" id="MathJax-Span-726"><span class="mn" id="MathJax-Span-727" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-728" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-729" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-730"><span class="mrow" id="MathJax-Span-731"><span class="mfrac" id="MathJax-Span-732"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-733"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-734" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-735"><span class="mrow" id="MathJax-Span-736"><span class="mn" id="MathJax-Span-737" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-738" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-739"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-740" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-741"><span class="mrow" id="MathJax-Span-742"><span class="mn" id="MathJax-Span-743" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-744" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-745" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-746"><span class="mrow" id="MathJax-Span-747"><span class="mfrac" id="MathJax-Span-748"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-749"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-750" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-751"><span class="mrow" id="MathJax-Span-752"><span class="mn" id="MathJax-Span-753" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-754" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-755"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-756" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-757"><span class="mrow" id="MathJax-Span-758"><span class="mn" id="MathJax-Span-759" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-760" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-761" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-762" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-763"><span class="mo" id="MathJax-Span-764" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-765"><span class="msubsup" id="MathJax-Span-766"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-767"><span class="mo" id="MathJax-Span-768" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-769"><span class="msubsup" id="MathJax-Span-770"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-771" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-772"><span class="mrow" id="MathJax-Span-773"><span class="mfrac" id="MathJax-Span-774"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-775"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-776" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-777"><span class="mrow" id="MathJax-Span-778"><span class="mn" id="MathJax-Span-779" style="font-size: 50%; font-family: MathJax_Main;">6</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-780"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-781" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-782"><span class="mrow" id="MathJax-Span-783"><span class="mn" id="MathJax-Span-784" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-785" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-786" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-787" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-788" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-789"><span class="mrow" id="MathJax-Span-790"><span class="mfrac" id="MathJax-Span-791"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-792"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-793" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-794"><span class="mrow" id="MathJax-Span-795"><span class="mn" id="MathJax-Span-796" style="font-size: 50%; font-family: MathJax_Main;">7</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-797"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-798" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-799"><span class="mrow" id="MathJax-Span-800"><span class="mn" id="MathJax-Span-801" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-802" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-803" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-804" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-805" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-806"><span class="mrow" id="MathJax-Span-807"><span class="mfrac" id="MathJax-Span-808"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-809"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-810" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-811"><span class="mrow" id="MathJax-Span-812"><span class="mn" id="MathJax-Span-813" style="font-size: 50%; font-family: MathJax_Main;">8</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-814"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-815" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-816"><span class="mrow" id="MathJax-Span-817"><span class="mn" id="MathJax-Span-818" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-819" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-820" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-821"><span class="mrow" id="MathJax-Span-822"><span class="mfrac" id="MathJax-Span-823"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-824"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-825" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-826"><span class="mrow" id="MathJax-Span-827"><span class="mn" id="MathJax-Span-828" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-829" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-830"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-831" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-832"><span class="mrow" id="MathJax-Span-833"><span class="mn" id="MathJax-Span-834" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-835" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-836" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-837" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-838"><span class="mo" id="MathJax-Span-839" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-840"><span class="msubsup" id="MathJax-Span-841"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-842" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-843"><span class="mrow" id="MathJax-Span-844"><span class="mfrac" id="MathJax-Span-845"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-846"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-847" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-848"><span class="mrow" id="MathJax-Span-849"><span class="mn" id="MathJax-Span-850" style="font-size: 50%; font-family: MathJax_Main;">10</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-851"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-852" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-853"><span class="mrow" id="MathJax-Span-854"><span class="mn" id="MathJax-Span-855" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-856" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-857" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-858" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-859" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-860"><span class="mrow" id="MathJax-Span-861"><span class="mfrac" id="MathJax-Span-862"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-863"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-864" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-865"><span class="mrow" id="MathJax-Span-866"><span class="mn" id="MathJax-Span-867" style="font-size: 50%; font-family: MathJax_Main;">11</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-868"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-869" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-870"><span class="mrow" id="MathJax-Span-871"><span class="mn" id="MathJax-Span-872" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-873" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-874" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-875" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-876" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-877"><span class="mrow" id="MathJax-Span-878"><span class="mfrac" id="MathJax-Span-879"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-880"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-881" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-882"><span class="mrow" id="MathJax-Span-883"><span class="mn" id="MathJax-Span-884" style="font-size: 50%; font-family: MathJax_Main;">9</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-885"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-886" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-887"><span class="mrow" id="MathJax-Span-888"><span class="mn" id="MathJax-Span-889" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-890" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-891" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 9.051em;"><span class="texatom" id="MathJax-Span-892"><span class="mrow" id="MathJax-Span-893"><span class="mfrac" id="MathJax-Span-894"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-895"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-896" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-897"><span class="mrow" id="MathJax-Span-898"><span class="mn" id="MathJax-Span-899" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-900" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-901"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-902" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-903"><span class="mrow" id="MathJax-Span-904"><span class="mn" id="MathJax-Span-905" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-906" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-907" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-908"><span class="mrow" id="MathJax-Span-909"><span class="mfrac" id="MathJax-Span-910"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-911"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-912" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-913"><span class="mrow" id="MathJax-Span-914"><span class="mn" id="MathJax-Span-915" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-916" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-917"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-918" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-919"><span class="mrow" id="MathJax-Span-920"><span class="mn" id="MathJax-Span-921" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-922" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-923" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 5.777em;"></span></span><span style="position: absolute; top: -6.366em; left: 52.324em;"><span class="texatom" id="MathJax-Span-924"><span class="mrow" id="MathJax-Span-925"><span class="msubsup" id="MathJax-Span-926"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.336em, 1000.3em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-927" style="font-size: 70.7%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.866em; left: 0.36em;"><span class="texatom" id="MathJax-Span-928"><span class="mrow" id="MathJax-Span-929"><span class="mn" id="MathJax-Span-930" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-931" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 6.253em;"></span></span><span style="position: absolute; clip: rect(3.217em, 1054.65em, 8.813em, -999.997em); top: -0.592em; left: 0em;"><span class="mrow" id="MathJax-Span-593-MathJax-Continue-1"><span class="mrow" id="MathJax-Span-595-MathJax-Continue-1"><span class="mo" id="MathJax-Span-932" style="font-family: MathJax_Main;">+</span><span class="mn" id="MathJax-Span-933" style="font-family: MathJax_Main; padding-left: 0.241em;">1</span></span><span class="mo" id="MathJax-Span-934" style="vertical-align: 2.86em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.598em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.604em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.068em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.473em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mrow" id="MathJax-Span-935" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-936" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎜<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-937"><span class="msubsup" id="MathJax-Span-938"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-939"><span class="mo" id="MathJax-Span-940" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-941"><span class="msubsup" id="MathJax-Span-942"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-943"><span class="mo" id="MathJax-Span-944" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-945"><span class="msubsup" id="MathJax-Span-946"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-947" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-948"><span class="mrow" id="MathJax-Span-949"><span class="mfrac" id="MathJax-Span-950"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-951"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-952" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-953"><span class="mrow" id="MathJax-Span-954"><span class="mn" id="MathJax-Span-955" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-956"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-957" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-958"><span class="mrow" id="MathJax-Span-959"><span class="mn" id="MathJax-Span-960" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-961" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-962" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-963" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-964" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-965"><span class="mrow" id="MathJax-Span-966"><span class="mfrac" id="MathJax-Span-967"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-968"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-969" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-970"><span class="mrow" id="MathJax-Span-971"><span class="mn" id="MathJax-Span-972" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-973"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-974" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-975"><span class="mrow" id="MathJax-Span-976"><span class="mn" id="MathJax-Span-977" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-978" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-979" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-980" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-981" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-982"><span class="mrow" id="MathJax-Span-983"><span class="mfrac" id="MathJax-Span-984"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-985"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-986" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-987"><span class="mrow" id="MathJax-Span-988"><span class="mn" id="MathJax-Span-989" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-990"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-991" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-992"><span class="mrow" id="MathJax-Span-993"><span class="mn" id="MathJax-Span-994" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-995" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-996" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-997"><span class="mrow" id="MathJax-Span-998"><span class="mfrac" id="MathJax-Span-999"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1000"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1001" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1002"><span class="mrow" id="MathJax-Span-1003"><span class="mn" id="MathJax-Span-1004" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1005" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1006"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1007" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1008"><span class="mrow" id="MathJax-Span-1009"><span class="mn" id="MathJax-Span-1010" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1011" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1012" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1013" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-1014"><span class="mo" id="MathJax-Span-1015" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1016"><span class="msubsup" id="MathJax-Span-1017"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1018" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1019"><span class="mrow" id="MathJax-Span-1020"><span class="mfrac" id="MathJax-Span-1021"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1022"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1023" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1024"><span class="mrow" id="MathJax-Span-1025"><span class="mn" id="MathJax-Span-1026" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1027"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1028" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1029"><span class="mrow" id="MathJax-Span-1030"><span class="mn" id="MathJax-Span-1031" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1032" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1033" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1034" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1035" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1036"><span class="mrow" id="MathJax-Span-1037"><span class="mfrac" id="MathJax-Span-1038"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1039"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1040" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1041"><span class="mrow" id="MathJax-Span-1042"><span class="mn" id="MathJax-Span-1043" style="font-size: 50%; font-family: MathJax_Main;">4</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1044"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1045" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1046"><span class="mrow" id="MathJax-Span-1047"><span class="mn" id="MathJax-Span-1048" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1049" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1050" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1051" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1052" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1053"><span class="mrow" id="MathJax-Span-1054"><span class="mfrac" id="MathJax-Span-1055"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1056"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1057" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1058"><span class="mrow" id="MathJax-Span-1059"><span class="mn" id="MathJax-Span-1060" style="font-size: 50%; font-family: MathJax_Main;">5</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1061"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1062" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1063"><span class="mrow" id="MathJax-Span-1064"><span class="mn" id="MathJax-Span-1065" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1066" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1067" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-1068"><span class="mrow" id="MathJax-Span-1069"><span class="mfrac" id="MathJax-Span-1070"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1071"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1072" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1073"><span class="mrow" id="MathJax-Span-1074"><span class="mn" id="MathJax-Span-1075" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1076" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1077"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1078" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1079"><span class="mrow" id="MathJax-Span-1080"><span class="mn" id="MathJax-Span-1081" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1082" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1083" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-1084"><span class="mrow" id="MathJax-Span-1085"><span class="mfrac" id="MathJax-Span-1086"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1087"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1088" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1089"><span class="mrow" id="MathJax-Span-1090"><span class="mn" id="MathJax-Span-1091" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1092" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1093"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1094" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1095"><span class="mrow" id="MathJax-Span-1096"><span class="mn" id="MathJax-Span-1097" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1098" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1099" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1100" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 24.646em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1023.16em, 6.313em, -999.997em); top: -4.58em; left: 0em;"><span class="mrow" id="MathJax-Span-1101"><span class="mo" id="MathJax-Span-1102" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-1103"><span class="msubsup" id="MathJax-Span-1104"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-1105"><span class="mo" id="MathJax-Span-1106" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1107"><span class="msubsup" id="MathJax-Span-1108"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1109" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1110"><span class="mrow" id="MathJax-Span-1111"><span class="mfrac" id="MathJax-Span-1112"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1113"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1114" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1115"><span class="mrow" id="MathJax-Span-1116"><span class="mn" id="MathJax-Span-1117" style="font-size: 50%; font-family: MathJax_Main;">6</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1118"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1119" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1120"><span class="mrow" id="MathJax-Span-1121"><span class="mn" id="MathJax-Span-1122" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1123" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1124" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1125" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1126" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1127"><span class="mrow" id="MathJax-Span-1128"><span class="mfrac" id="MathJax-Span-1129"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1130"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1131" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1132"><span class="mrow" id="MathJax-Span-1133"><span class="mn" id="MathJax-Span-1134" style="font-size: 50%; font-family: MathJax_Main;">7</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1135"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1136" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1137"><span class="mrow" id="MathJax-Span-1138"><span class="mn" id="MathJax-Span-1139" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1140" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1141" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1142" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1143" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1144"><span class="mrow" id="MathJax-Span-1145"><span class="mfrac" id="MathJax-Span-1146"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1147"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1148" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1149"><span class="mrow" id="MathJax-Span-1150"><span class="mn" id="MathJax-Span-1151" style="font-size: 50%; font-family: MathJax_Main;">8</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1152"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1153" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1154"><span class="mrow" id="MathJax-Span-1155"><span class="mn" id="MathJax-Span-1156" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1157" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1158" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-1159"><span class="mrow" id="MathJax-Span-1160"><span class="mfrac" id="MathJax-Span-1161"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1162"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1163" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1164"><span class="mrow" id="MathJax-Span-1165"><span class="mn" id="MathJax-Span-1166" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1167" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1168"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1169" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1170"><span class="mrow" id="MathJax-Span-1171"><span class="mn" id="MathJax-Span-1172" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1173" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1174" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1175" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-1176"><span class="mo" id="MathJax-Span-1177" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1178"><span class="msubsup" id="MathJax-Span-1179"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1180" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1181"><span class="mrow" id="MathJax-Span-1182"><span class="mfrac" id="MathJax-Span-1183"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-1184"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1185" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1186"><span class="mrow" id="MathJax-Span-1187"><span class="mn" id="MathJax-Span-1188" style="font-size: 50%; font-family: MathJax_Main;">10</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1189"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1190" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1191"><span class="mrow" id="MathJax-Span-1192"><span class="mn" id="MathJax-Span-1193" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1194" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1195" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1196" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1197" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1198"><span class="mrow" id="MathJax-Span-1199"><span class="mfrac" id="MathJax-Span-1200"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.414em;"><span class="msubsup" id="MathJax-Span-1201"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1202" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1203"><span class="mrow" id="MathJax-Span-1204"><span class="mn" id="MathJax-Span-1205" style="font-size: 50%; font-family: MathJax_Main;">11</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1206"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1207" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1208"><span class="mrow" id="MathJax-Span-1209"><span class="mn" id="MathJax-Span-1210" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1211" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1212" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1213" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1214" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1215"><span class="mrow" id="MathJax-Span-1216"><span class="mfrac" id="MathJax-Span-1217"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1218"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1219" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1220"><span class="mrow" id="MathJax-Span-1221"><span class="mn" id="MathJax-Span-1222" style="font-size: 50%; font-family: MathJax_Main;">9</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1223"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1224" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1225"><span class="mrow" id="MathJax-Span-1226"><span class="mn" id="MathJax-Span-1227" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1228" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1229" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 9.051em;"><span class="texatom" id="MathJax-Span-1230"><span class="mrow" id="MathJax-Span-1231"><span class="mfrac" id="MathJax-Span-1232"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1233"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1234" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1235"><span class="mrow" id="MathJax-Span-1236"><span class="mn" id="MathJax-Span-1237" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1238" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1239"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1240" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1241"><span class="mrow" id="MathJax-Span-1242"><span class="mn" id="MathJax-Span-1243" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1244" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1245" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span><span style="position: absolute; top: -5.771em; left: 23.455em;"><span class="texatom" id="MathJax-Span-1246"><span class="mrow" id="MathJax-Span-1247"><span class="mfrac" id="MathJax-Span-1248"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1249"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1250" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1251"><span class="mrow" id="MathJax-Span-1252"><span class="mn" id="MathJax-Span-1253" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1254" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1255"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1256" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1257"><span class="mrow" id="MathJax-Span-1258"><span class="mn" id="MathJax-Span-1259" style="font-size: 50%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1260" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1261" style="vertical-align: 2.622em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: 0.182em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.723em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -1.307em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="font-family: MathJax_Size4; position: absolute; top: -0.83em; left: 0em;">⎟<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mrow" id="MathJax-Span-1262" style="padding-left: 0.182em;"></span><span style="display: inline-block; width: 0px; height: 6.253em;"></span></span><span style="position: absolute; clip: rect(2.384em, 1032.44em, 6.313em, -999.997em); top: 5.955em; left: 0em;"><span class="mrow" id="MathJax-Span-1262-MathJax-Continue-1" style=""><span class="mo" id="MathJax-Span-1263" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎛<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎝<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-1264"><span class="msubsup" id="MathJax-Span-1265"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-1266"><span class="mo" id="MathJax-Span-1267" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1268"><span class="msubsup" id="MathJax-Span-1269"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1270" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1271"><span class="mrow" id="MathJax-Span-1272"><span class="mfrac" id="MathJax-Span-1273"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1274"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1275" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1276"><span class="mrow" id="MathJax-Span-1277"><span class="mn" id="MathJax-Span-1278" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1279"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1280" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1281"><span class="mrow" id="MathJax-Span-1282"><span class="mn" id="MathJax-Span-1283" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1284" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1285" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1286" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1287" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1288"><span class="mrow" id="MathJax-Span-1289"><span class="mfrac" id="MathJax-Span-1290"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1291"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1292" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1293"><span class="mrow" id="MathJax-Span-1294"><span class="mn" id="MathJax-Span-1295" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1296"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1297" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1298"><span class="mrow" id="MathJax-Span-1299"><span class="mn" id="MathJax-Span-1300" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1301" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1302" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1303" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1304" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1305"><span class="mrow" id="MathJax-Span-1306"><span class="mfrac" id="MathJax-Span-1307"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1308"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1309" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1310"><span class="mrow" id="MathJax-Span-1311"><span class="mn" id="MathJax-Span-1312" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1313"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1314" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1315"><span class="mrow" id="MathJax-Span-1316"><span class="mn" id="MathJax-Span-1317" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1318" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1319" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-1320"><span class="mrow" id="MathJax-Span-1321"><span class="mfrac" id="MathJax-Span-1322"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1323"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1324" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1325"><span class="mrow" id="MathJax-Span-1326"><span class="mn" id="MathJax-Span-1327" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1328" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1329"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1330" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1331"><span class="mrow" id="MathJax-Span-1332"><span class="mn" id="MathJax-Span-1333" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1334" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1335" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1336" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 10.241em; height: 0px;"><span style="position: absolute; clip: rect(2.384em, 1008.81em, 5.122em, -999.997em); top: -3.985em; left: 0em;"><span class="mrow" id="MathJax-Span-1337"><span class="mo" id="MathJax-Span-1338" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1339"><span class="msubsup" id="MathJax-Span-1340"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1341" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1342"><span class="mrow" id="MathJax-Span-1343"><span class="mfrac" id="MathJax-Span-1344"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1345"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1346" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1347"><span class="mrow" id="MathJax-Span-1348"><span class="mn" id="MathJax-Span-1349" style="font-size: 50%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1350"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1351" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1352"><span class="mrow" id="MathJax-Span-1353"><span class="mn" id="MathJax-Span-1354" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1355" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1356" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1357" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1358" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1359"><span class="mrow" id="MathJax-Span-1360"><span class="mfrac" id="MathJax-Span-1361"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1362"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1363" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1364"><span class="mrow" id="MathJax-Span-1365"><span class="mn" id="MathJax-Span-1366" style="font-size: 50%; font-family: MathJax_Main;">4</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1367"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1368" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1369"><span class="mrow" id="MathJax-Span-1370"><span class="mn" id="MathJax-Span-1371" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1372" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1373" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1374" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1375" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.58em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1376"><span class="mrow" id="MathJax-Span-1377"><span class="mfrac" id="MathJax-Span-1378"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1379"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1380" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1381"><span class="mrow" id="MathJax-Span-1382"><span class="mn" id="MathJax-Span-1383" style="font-size: 50%; font-family: MathJax_Main;">5</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1384"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1385" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1386"><span class="mrow" id="MathJax-Span-1387"><span class="mn" id="MathJax-Span-1388" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1389" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1390" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -5.176em; left: 8.991em;"><span class="texatom" id="MathJax-Span-1391"><span class="mrow" id="MathJax-Span-1392"><span class="mfrac" id="MathJax-Span-1393"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1394"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1395" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1396"><span class="mrow" id="MathJax-Span-1397"><span class="mn" id="MathJax-Span-1398" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1399" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1400"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1401" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1402"><span class="mrow" id="MathJax-Span-1403"><span class="mn" id="MathJax-Span-1404" style="font-size: 50%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1405" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1406" style="vertical-align: 2.027em;"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; font-family: MathJax_Size4; top: -2.854em; left: 0em;">⎞<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; font-family: MathJax_Size4; top: -1.009em; left: 0em;">⎠<span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mrow" id="MathJax-Span-1407" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-1408" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1409"><span class="msubsup" id="MathJax-Span-1410"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1411" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1412"><span class="mrow" id="MathJax-Span-1413"><span class="mfrac" id="MathJax-Span-1414"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1415"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1416" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1417"><span class="mrow" id="MathJax-Span-1418"><span class="mn" id="MathJax-Span-1419" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1420"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1421" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1422"><span class="mrow" id="MathJax-Span-1423"><span class="mn" id="MathJax-Span-1424" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1425" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1426" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1427" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1428" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1429"><span class="mrow" id="MathJax-Span-1430"><span class="mfrac" id="MathJax-Span-1431"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1432"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1433" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1434"><span class="mrow" id="MathJax-Span-1435"><span class="mn" id="MathJax-Span-1436" style="font-size: 50%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1437"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1438" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1439"><span class="mrow" id="MathJax-Span-1440"><span class="mn" id="MathJax-Span-1441" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1442" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1443" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="msubsup" id="MathJax-Span-1444" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.67em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1445" style="font-family: MathJax_Math-italic;">e</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -4.64em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1446"><span class="mrow" id="MathJax-Span-1447"><span class="mfrac" id="MathJax-Span-1448"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.574em, 1000.54em, 4.348em, -999.997em); top: -4.521em; left: 50%; margin-left: -0.235em;"><span class="msubsup" id="MathJax-Span-1449"><span style="display: inline-block; position: relative; width: 0.539em; height: 0px;"><span style="position: absolute; clip: rect(3.574em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1450" style="font-size: 50%; font-family: MathJax_Math-italic;">v</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1451"><span class="mrow" id="MathJax-Span-1452"><span class="mn" id="MathJax-Span-1453" style="font-size: 50%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.455em, 1000.78em, 4.348em, -999.997em); top: -3.628em; left: 50%; margin-left: -0.354em;"><span class="msubsup" id="MathJax-Span-1454"><span style="display: inline-block; position: relative; width: 0.777em; height: 0px;"><span style="position: absolute; clip: rect(3.455em, 1000.24em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1455" style="font-size: 50%; font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.241em;"><span class="texatom" id="MathJax-Span-1456"><span class="mrow" id="MathJax-Span-1457"><span class="mn" id="MathJax-Span-1458" style="font-size: 50%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1459" style="font-size: 50%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1000.9em, 1.253em, -999.997em); top: -1.188em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 0.896em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1460" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span style="display: inline-block; width: 0px; height: 4.586em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 15.955em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1055.24em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 55.241em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1461" style="font-family: MathJax_Main;">.</span></span><span style="display: inline-block; width: 0px; height: 2.205em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -17.996em; border-left: 0px solid; width: 0px; height: 30.218em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><msub><mi>s</mi><mi>j</mi></msub><mo>=</mo><mfrac><mrow><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>4</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>5</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>6</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>7</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>8</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>10</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>11</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>9</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></msup><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>4</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>5</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mrow><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>4</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>5</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>6</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>7</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>8</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>10</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>11</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>9</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></msup><mo>+</mo><mn>1</mn></mrow><mo>)</mo></mrow><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>4</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>5</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>6</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>7</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>8</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>10</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>11</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>9</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow><mo>(</mo><mrow><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>4</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>5</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow><mrow><mo>(</mo><mrow><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup><mo>+</mo><msup><mi>e</mi><mrow class="MJX-TeXAtom-ORD"><mfrac><msub><mi>v</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac></mrow></msup></mrow><mo>)</mo></mrow></mrow></mfrac><mo>.</mo></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-1"> s_j = \frac{\left(\left(\left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} + \left(e^{\frac{v_{3}}{\theta_{2 1}}} + e^{\frac{v_{4}}{\theta_{2 1}}} + e^{\frac{v_{5}}{\theta_{2 1}}}\right)^{\frac{\theta_{2 1}}{\theta_{1 0}}}\right)^{\frac{\theta_{1 0}}{\theta_{0 0}}} + \left(\left(e^{\frac{v_{6}}{\theta_{2 2}}} + e^{\frac{v_{7}}{\theta_{2 2}}} + e^{\frac{v_{8}}{\theta_{2 2}}}\right)^{\frac{\theta_{2 2}}{\theta_{1 1}}} + \left(e^{\frac{v_{10}}{\theta_{2 3}}} + e^{\frac{v_{11}}{\theta_{2 3}}} + e^{\frac{v_{9}}{\theta_{2 3}}}\right)^{\frac{\theta_{2 3}}{\theta_{1 1}}}\right)^{\frac{\theta_{1 1}}{\theta_{0 0}}}\right)^{\theta_{0 0}} \left(\left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} + \left(e^{\frac{v_{3}}{\theta_{2 1}}} + e^{\frac{v_{4}}{\theta_{2 1}}} + e^{\frac{v_{5}}{\theta_{2 1}}}\right)^{\frac{\theta_{2 1}}{\theta_{1 0}}}\right)^{\frac{\theta_{1 0}}{\theta_{0 0}}} \left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} e^{\frac{v_{0}}{\theta_{2 0}}}}{\left(\left(\left(\left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} + \left(e^{\frac{v_{3}}{\theta_{2 1}}} + e^{\frac{v_{4}}{\theta_{2 1}}} + e^{\frac{v_{5}}{\theta_{2 1}}}\right)^{\frac{\theta_{2 1}}{\theta_{1 0}}}\right)^{\frac{\theta_{1 0}}{\theta_{0 0}}} + \left(\left(e^{\frac{v_{6}}{\theta_{2 2}}} + e^{\frac{v_{7}}{\theta_{2 2}}} + e^{\frac{v_{8}}{\theta_{2 2}}}\right)^{\frac{\theta_{2 2}}{\theta_{1 1}}} + \left(e^{\frac{v_{10}}{\theta_{2 3}}} + e^{\frac{v_{11}}{\theta_{2 3}}} + e^{\frac{v_{9}}{\theta_{2 3}}}\right)^{\frac{\theta_{2 3}}{\theta_{1 1}}}\right)^{\frac{\theta_{1 1}}{\theta_{0 0}}}\right)^{\theta_{0 0}} + 1\right) \left(\left(\left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} + \left(e^{\frac{v_{3}}{\theta_{2 1}}} + e^{\frac{v_{4}}{\theta_{2 1}}} + e^{\frac{v_{5}}{\theta_{2 1}}}\right)^{\frac{\theta_{2 1}}{\theta_{1 0}}}\right)^{\frac{\theta_{1 0}}{\theta_{0 0}}} + \left(\left(e^{\frac{v_{6}}{\theta_{2 2}}} + e^{\frac{v_{7}}{\theta_{2 2}}} + e^{\frac{v_{8}}{\theta_{2 2}}}\right)^{\frac{\theta_{2 2}}{\theta_{1 1}}} + \left(e^{\frac{v_{10}}{\theta_{2 3}}} + e^{\frac{v_{11}}{\theta_{2 3}}} + e^{\frac{v_{9}}{\theta_{2 3}}}\right)^{\frac{\theta_{2 3}}{\theta_{1 1}}}\right)^{\frac{\theta_{1 1}}{\theta_{0 0}}}\right) \left(\left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)^{\frac{\theta_{2 0}}{\theta_{1 0}}} + \left(e^{\frac{v_{3}}{\theta_{2 1}}} + e^{\frac{v_{4}}{\theta_{2 1}}} + e^{\frac{v_{5}}{\theta_{2 1}}}\right)^{\frac{\theta_{2 1}}{\theta_{1 0}}}\right) \left(e^{\frac{v_{0}}{\theta_{2 0}}} + e^{\frac{v_{1}}{\theta_{2 0}}} + e^{\frac{v_{2}}{\theta_{2 0}}}\right)}. </script></p>

</div>
</div>
</div>
<div class="cell border-box-sizing text_cell rendered"><div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<h2 id="Self-elasticity">Self-elasticity<a class="anchor-link" href="http://localhost:8888/nbconvert/html/mdm_model_formula_check.ipynb?download=false#Self-elasticity">¶</a></h2><p>The self-elasticity of the node is
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-2-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mi&gt;l&lt;/mi&gt;&lt;mi&gt;f&lt;/mi&gt;&lt;mi mathvariant=&quot;normal&quot;&gt;&amp;#x005F;&lt;/mi&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mi&gt;l&lt;/mi&gt;&lt;mi&gt;a&lt;/mi&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mo&gt;=&lt;/mo&gt;&lt;mi&gt;&amp;#x03B2;&lt;/mi&gt;&lt;msub&gt;&lt;mi&gt;p&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mo&gt;,&lt;/mo&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-1462" style="width: 41.967em; display: inline-block;"><span style="display: inline-block; position: relative; width: 34.943em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(1.789em, 1034.88em, 4.527em, -999.997em); top: -3.39em; left: 0em;"><span class="mrow" id="MathJax-Span-1463"><span class="mi" id="MathJax-Span-1464" style="font-family: MathJax_Math-italic;">s</span><span class="mi" id="MathJax-Span-1465" style="font-family: MathJax_Math-italic;">e</span><span class="mi" id="MathJax-Span-1466" style="font-family: MathJax_Math-italic;">l</span><span class="mi" id="MathJax-Span-1467" style="font-family: MathJax_Math-italic;">f<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.063em;"></span></span><span class="mi" id="MathJax-Span-1468" style="font-family: MathJax_Main;">_</span><span class="mi" id="MathJax-Span-1469" style="font-family: MathJax_Math-italic;">e</span><span class="mi" id="MathJax-Span-1470" style="font-family: MathJax_Math-italic;">l</span><span class="mi" id="MathJax-Span-1471" style="font-family: MathJax_Math-italic;">a</span><span class="mi" id="MathJax-Span-1472" style="font-family: MathJax_Math-italic;">s</span><span class="mo" id="MathJax-Span-1473" style="font-family: MathJax_Main; padding-left: 0.301em;">=</span><span class="mi" id="MathJax-Span-1474" style="font-family: MathJax_Math-italic; padding-left: 0.301em;">β<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.003em;"></span></span><span class="msubsup" id="MathJax-Span-1475"><span style="display: inline-block; position: relative; width: 0.955em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.48em, 4.348em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1476" style="font-family: MathJax_Math-italic;">p</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1477"><span class="mrow" id="MathJax-Span-1478"><span class="mn" id="MathJax-Span-1479" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-1480" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-1481" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1482"><span class="mo" id="MathJax-Span-1483" style="font-family: MathJax_Main;">−</span><span class="msubsup" id="MathJax-Span-1484"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1485" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1486"><span class="mrow" id="MathJax-Span-1487"><span class="mn" id="MathJax-Span-1488" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1489" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1490" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.372em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.158em, 1000.42em, 4.17em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.235em;"><span class="mn" id="MathJax-Span-1491" style="font-family: MathJax_Main;">1</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1001.25em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.652em;"><span class="msubsup" id="MathJax-Span-1492"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1493" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1494"><span class="mrow" id="MathJax-Span-1495"><span class="mn" id="MathJax-Span-1496" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1497" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.37em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.372em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1498" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1499" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1500"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1501" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1502"><span class="mrow" id="MathJax-Span-1503"><span class="mn" id="MathJax-Span-1504" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1505"><span class="msubsup" id="MathJax-Span-1506"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1507" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1508"><span class="mrow" id="MathJax-Span-1509"><span class="mn" id="MathJax-Span-1510" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1511"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1512" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1513"><span class="mrow" id="MathJax-Span-1514"><span class="mn" id="MathJax-Span-1515" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1516" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1517" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1518" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1519"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1520" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1521"><span class="mrow" id="MathJax-Span-1522"><span class="mn" id="MathJax-Span-1523" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1524"><span class="msubsup" id="MathJax-Span-1525"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1526" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1527"><span class="mrow" id="MathJax-Span-1528"><span class="mn" id="MathJax-Span-1529" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1530"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1531" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1532"><span class="mrow" id="MathJax-Span-1533"><span class="mn" id="MathJax-Span-1534" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1535" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1536" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1537" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1538"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1539" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1540"><span class="mrow" id="MathJax-Span-1541"><span class="mn" id="MathJax-Span-1542" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1543"><span class="msubsup" id="MathJax-Span-1544"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1545" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1546"><span class="mrow" id="MathJax-Span-1547"><span class="mn" id="MathJax-Span-1548" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1549"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1550" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1551"><span class="mrow" id="MathJax-Span-1552"><span class="mn" id="MathJax-Span-1553" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1554" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1555" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1556" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1557"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1558" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1559"><span class="mrow" id="MathJax-Span-1560"><span class="mn" id="MathJax-Span-1561" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1562"><span class="msubsup" id="MathJax-Span-1563"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1564" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1565"><span class="mrow" id="MathJax-Span-1566"><span class="mn" id="MathJax-Span-1567" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1568"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1569" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1570"><span class="mrow" id="MathJax-Span-1571"><span class="mn" id="MathJax-Span-1572" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1573" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1574" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1575" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.015em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1576"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1577" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1578"><span class="mrow" id="MathJax-Span-1579"><span class="mn" id="MathJax-Span-1580" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1581"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1582" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1583"><span class="mrow" id="MathJax-Span-1584"><span class="mn" id="MathJax-Span-1585" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.01em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.015em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1586" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1587" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1588"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1589" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1590"><span class="mrow" id="MathJax-Span-1591"><span class="mn" id="MathJax-Span-1592" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1593"><span class="msubsup" id="MathJax-Span-1594"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1595" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1596"><span class="mrow" id="MathJax-Span-1597"><span class="mn" id="MathJax-Span-1598" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1599"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1600" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1601"><span class="mrow" id="MathJax-Span-1602"><span class="mn" id="MathJax-Span-1603" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1604" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1605" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span class="mo" id="MathJax-Span-1606" style="font-family: MathJax_Main; padding-left: 0.182em;">,</span></span><span style="display: inline-block; width: 0px; height: 3.396em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -1.211em; border-left: 0px solid; width: 0px; height: 3.004em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>s</mi><mi>e</mi><mi>l</mi><mi>f</mi><mi mathvariant="normal">_</mi><mi>e</mi><mi>l</mi><mi>a</mi><mi>s</mi><mo>=</mo><mi>β</mi><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><mrow><mo>(</mo><mrow><mo>−</mo><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mo>+</mo><mfrac><mn>1</mn><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac></mrow><mo>)</mo></mrow><mo>,</mo></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-2"> self\_elas = \beta p_{0} \left(- s_{3} + \frac{1}{\theta_{2 0}} - \frac{s_{3}}{s_{2} \theta_{2 0}} + \frac{s_{3}}{s_{2} \theta_{1 0}} - \frac{s_{3}}{s_{1} \theta_{1 0}} + \frac{s_{3}}{s_{1} \theta_{0 0}} + \frac{s_{3}}{s_{0}} - \frac{s_{3}}{s_{0} \theta_{0 0}}\right), </script> and the cross-elasticity of the node with another node within the same group in the second to last layer is
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-3-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mi&gt;c&lt;/mi&gt;&lt;mi&gt;r&lt;/mi&gt;&lt;mi&gt;o&lt;/mi&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mi mathvariant=&quot;normal&quot;&gt;&amp;#x005F;&lt;/mi&gt;&lt;mi&gt;e&lt;/mi&gt;&lt;mi&gt;l&lt;/mi&gt;&lt;mi&gt;a&lt;/mi&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mo&gt;=&lt;/mo&gt;&lt;mi&gt;&amp;#x03B2;&lt;/mi&gt;&lt;msub&gt;&lt;mi&gt;p&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;mo&gt;.&lt;/mo&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-1607" style="width: 39.17em; display: inline-block;"><span style="display: inline-block; position: relative; width: 32.622em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(1.789em, 1032.56em, 4.527em, -999.997em); top: -3.39em; left: 0em;"><span class="mrow" id="MathJax-Span-1608"><span class="mi" id="MathJax-Span-1609" style="font-family: MathJax_Math-italic;">c</span><span class="mi" id="MathJax-Span-1610" style="font-family: MathJax_Math-italic;">r</span><span class="mi" id="MathJax-Span-1611" style="font-family: MathJax_Math-italic;">o</span><span class="mi" id="MathJax-Span-1612" style="font-family: MathJax_Math-italic;">s</span><span class="mi" id="MathJax-Span-1613" style="font-family: MathJax_Math-italic;">s</span><span class="mi" id="MathJax-Span-1614" style="font-family: MathJax_Main;">_</span><span class="mi" id="MathJax-Span-1615" style="font-family: MathJax_Math-italic;">e</span><span class="mi" id="MathJax-Span-1616" style="font-family: MathJax_Math-italic;">l</span><span class="mi" id="MathJax-Span-1617" style="font-family: MathJax_Math-italic;">a</span><span class="mi" id="MathJax-Span-1618" style="font-family: MathJax_Math-italic;">s</span><span class="mo" id="MathJax-Span-1619" style="font-family: MathJax_Main; padding-left: 0.301em;">=</span><span class="mi" id="MathJax-Span-1620" style="font-family: MathJax_Math-italic; padding-left: 0.301em;">β<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.003em;"></span></span><span class="msubsup" id="MathJax-Span-1621"><span style="display: inline-block; position: relative; width: 0.955em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.48em, 4.348em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1622" style="font-family: MathJax_Math-italic;">p</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1623"><span class="mrow" id="MathJax-Span-1624"><span class="mn" id="MathJax-Span-1625" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-1626" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-1627" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-1628"><span class="mo" id="MathJax-Span-1629" style="font-family: MathJax_Main;">−</span><span class="msubsup" id="MathJax-Span-1630"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1631" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1632"><span class="mrow" id="MathJax-Span-1633"><span class="mn" id="MathJax-Span-1634" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-1635" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1636" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1637"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1638" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1639"><span class="mrow" id="MathJax-Span-1640"><span class="mn" id="MathJax-Span-1641" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1642"><span class="msubsup" id="MathJax-Span-1643"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1644" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1645"><span class="mrow" id="MathJax-Span-1646"><span class="mn" id="MathJax-Span-1647" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1648"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1649" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1650"><span class="mrow" id="MathJax-Span-1651"><span class="mn" id="MathJax-Span-1652" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-1653" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1654" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1655" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1656"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1657" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1658"><span class="mrow" id="MathJax-Span-1659"><span class="mn" id="MathJax-Span-1660" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1661"><span class="msubsup" id="MathJax-Span-1662"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1663" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1664"><span class="mrow" id="MathJax-Span-1665"><span class="mn" id="MathJax-Span-1666" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1667"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1668" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1669"><span class="mrow" id="MathJax-Span-1670"><span class="mn" id="MathJax-Span-1671" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1672" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1673" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1674" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1675"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1676" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1677"><span class="mrow" id="MathJax-Span-1678"><span class="mn" id="MathJax-Span-1679" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1680"><span class="msubsup" id="MathJax-Span-1681"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1682" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1683"><span class="mrow" id="MathJax-Span-1684"><span class="mn" id="MathJax-Span-1685" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1686"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1687" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1688"><span class="mrow" id="MathJax-Span-1689"><span class="mn" id="MathJax-Span-1690" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-1691" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1692" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1693" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1694"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1695" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1696"><span class="mrow" id="MathJax-Span-1697"><span class="mn" id="MathJax-Span-1698" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1699"><span class="msubsup" id="MathJax-Span-1700"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1701" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1702"><span class="mrow" id="MathJax-Span-1703"><span class="mn" id="MathJax-Span-1704" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1705"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1706" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1707"><span class="mrow" id="MathJax-Span-1708"><span class="mn" id="MathJax-Span-1709" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1710" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1711" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-1712" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.015em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1713"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1714" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1715"><span class="mrow" id="MathJax-Span-1716"><span class="mn" id="MathJax-Span-1717" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1718"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1719" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1720"><span class="mrow" id="MathJax-Span-1721"><span class="mn" id="MathJax-Span-1722" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.01em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.015em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-1723" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-1724" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-1725"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1726" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1727"><span class="mrow" id="MathJax-Span-1728"><span class="mn" id="MathJax-Span-1729" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-1730"><span class="msubsup" id="MathJax-Span-1731"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1732" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1733"><span class="mrow" id="MathJax-Span-1734"><span class="mn" id="MathJax-Span-1735" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-1736"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-1737" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-1738"><span class="mrow" id="MathJax-Span-1739"><span class="mn" id="MathJax-Span-1740" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-1741" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-1742" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span><span class="mo" id="MathJax-Span-1743" style="font-family: MathJax_Main; padding-left: 0.182em;">.</span></span><span style="display: inline-block; width: 0px; height: 3.396em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -1.211em; border-left: 0px solid; width: 0px; height: 3.004em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>c</mi><mi>r</mi><mi>o</mi><mi>s</mi><mi>s</mi><mi mathvariant="normal">_</mi><mi>e</mi><mi>l</mi><mi>a</mi><mi>s</mi><mo>=</mo><mi>β</mi><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><mrow><mo>(</mo><mrow><mo>−</mo><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac></mrow><mo>)</mo></mrow><mo>.</mo></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-3"> cross\_elas = \beta p_{1} \left(- s_{3} - \frac{s_{3}}{s_{2} \theta_{2 0}} + \frac{s_{3}}{s_{2} \theta_{1 0}} - \frac{s_{3}}{s_{1} \theta_{1 0}} + \frac{s_{3}}{s_{1} \theta_{0 0}} + \frac{s_{3}}{s_{0}} - \frac{s_{3}}{s_{0} \theta_{0 0}}\right). </script></p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[4]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1"># Basic elements for calculation</span>
<span class="n">total_end_nodes</span> <span class="o">=</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="p">(</span><span class="n">depth_of_tree</span> <span class="o">-</span> <span class="mi">1</span><span class="p">)</span> <span class="o">*</span> <span class="n">number_end_group</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">total_end_nodes</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'p_</span><span class="si">%d</span><span class="s1"> = symbols("p_</span><span class="si">%d</span><span class="s1">")'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="n">total_theta</span> <span class="o">=</span> <span class="mi">0</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="n">total_theta</span> <span class="o">+=</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="n">i</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">total_end_nodes</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'v_</span><span class="si">%d</span><span class="s1"> = symbols("v_</span><span class="si">%d</span><span class="s1">")'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="k">for</span> <span class="n">j</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="n">i</span><span class="p">):</span>
        <span class="n">exec</span><span class="p">(</span><span class="s1">'theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> = symbols("theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">")'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">))</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="o">+</span><span class="mi">1</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_</span><span class="si">%d</span><span class="s1"> = symbols("s_</span><span class="si">%d</span><span class="s1">")'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[5]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1"># Intermediate sums in s</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">total_end_nodes</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> = exp(v_</span><span class="si">%d</span><span class="s1"> / theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">))</span>
    

<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">total_end_nodes</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">):</span>
    <span class="n">st</span> <span class="o">=</span> <span class="s1">'sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> ='</span> <span class="o">%</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">)</span>
    <span class="k">for</span> <span class="n">j</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">i</span><span class="o">*</span><span class="n">number_end_group</span><span class="p">,</span> <span class="n">i</span><span class="o">*</span><span class="n">number_end_group</span> <span class="o">+</span> <span class="n">number_end_group</span><span class="p">):</span>
        <span class="n">st</span> <span class="o">=</span> <span class="n">st</span> <span class="o">+</span> <span class="s1">' exp(v_</span><span class="si">%d</span><span class="s1"> / theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">) +'</span> <span class="o">%</span> <span class="p">(</span><span class="n">j</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">)</span>
    <span class="n">exec</span><span class="p">(</span><span class="n">st</span><span class="p">[:</span><span class="o">-</span><span class="mi">1</span><span class="p">])</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> = sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> ** (theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> / theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span>
                                                                       <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">2</span><span class="p">,</span> <span class="n">i</span><span class="o">//</span><span class="n">number_splits_non_end_node</span><span class="p">))</span>

<span class="k">for</span> <span class="n">k</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">2</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="o">-</span><span class="mi">1</span><span class="p">):</span>
    <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="n">k</span><span class="p">):</span>
        <span class="n">st</span> <span class="o">=</span> <span class="s1">'sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> ='</span> <span class="o">%</span> <span class="p">(</span><span class="n">k</span><span class="p">,</span> <span class="n">i</span><span class="p">)</span>
        <span class="k">for</span> <span class="n">j</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">i</span><span class="o">*</span><span class="n">number_splits_non_end_node</span><span class="p">,</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">)</span><span class="o">*</span><span class="n">number_splits_non_end_node</span><span class="p">):</span>
            <span class="n">st</span> <span class="o">=</span> <span class="n">st</span> <span class="o">+</span> <span class="s1">' sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> +'</span> <span class="o">%</span> <span class="p">(</span><span class="n">k</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">j</span><span class="p">)</span>
        <span class="n">exec</span><span class="p">(</span><span class="n">st</span><span class="p">[:</span><span class="o">-</span><span class="mi">1</span><span class="p">])</span>
        <span class="n">exec</span><span class="p">(</span><span class="s1">'sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> = sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> ** (theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> / theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">k</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span>
                                                                       <span class="n">k</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="o">//</span><span class="n">number_splits_non_end_node</span><span class="p">))</span>

<span class="n">st</span> <span class="o">=</span> <span class="s1">'sum_0_0 ='</span>
<span class="k">for</span> <span class="n">j</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">number_splits_non_end_node</span><span class="p">):</span>
    <span class="n">st</span> <span class="o">=</span> <span class="n">st</span> <span class="o">+</span> <span class="s1">' sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> +'</span> <span class="o">%</span> <span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">j</span><span class="p">)</span>
<span class="n">exec</span><span class="p">(</span><span class="n">st</span><span class="p">[:</span><span class="o">-</span><span class="mi">1</span><span class="p">])</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'sum_exp_0_0 = sum_0_0 ** (theta_0_0)'</span><span class="p">)</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[6]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1"># Calculate final expression for s</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'s = sum_exp_</span><span class="si">%d</span><span class="s1">_0 / (1 + sum_0_0 ** theta_0_0)'</span> <span class="o">%</span> <span class="n">depth_of_tree</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s = s / sum_</span><span class="si">%d</span><span class="s1">_0 * sum_exp_</span><span class="si">%d</span><span class="s1">_0'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[7]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">s</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt output_prompt">Out[7]:</div>




<div class="output_latex output_subarea output_execute_result">
</div>

</div>

</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[8]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#### Deriving the own-elasticity formular</span>
<span class="n">s_diff</span> <span class="o">=</span> <span class="n">diff</span><span class="p">(</span><span class="n">s</span><span class="p">,</span> <span class="n">v_0</span><span class="p">)</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[9]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">s_simp</span> <span class="o">=</span> <span class="n">s_diff</span><span class="o">.</span><span class="n">subs</span><span class="p">(</span><span class="n">sum_exp_0_0</span> <span class="o">/</span> <span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">sum_exp_0_0</span><span class="p">),</span> <span class="n">s_0</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_simp = s_simp.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_0 / sum_</span><span class="si">%d</span><span class="s1">_0, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="n">s_simp</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt output_prompt">Out[9]:</div>




<div class="output_latex output_subarea output_execute_result">
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-5-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msubsup&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msubsup&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-3200" style="width: 31.074em; display: inline-block;"><span style="display: inline-block; position: relative; width: 25.896em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(0.42em, 1025.9em, 3.217em, -999.997em); top: -2.199em; left: 0em;"><span class="mrow" id="MathJax-Span-3201"><span class="mo" id="MathJax-Span-3202" style="font-family: MathJax_Main;">−</span><span class="msubsup" id="MathJax-Span-3203"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3204" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3205"><span class="mrow" id="MathJax-Span-3206"><span class="mn" id="MathJax-Span-3207" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3208"><span class="mrow" id="MathJax-Span-3209"><span class="mn" id="MathJax-Span-3210" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-3211" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3212" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.372em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3213"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3214" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3215"><span class="mrow" id="MathJax-Span-3216"><span class="mn" id="MathJax-Span-3217" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1001.25em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.652em;"><span class="msubsup" id="MathJax-Span-3218"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3219" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3220"><span class="mrow" id="MathJax-Span-3221"><span class="mn" id="MathJax-Span-3222" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-3223" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.37em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.372em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3224" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3225" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3226"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3227" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3228"><span class="mrow" id="MathJax-Span-3229"><span class="mn" id="MathJax-Span-3230" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3231"><span class="mrow" id="MathJax-Span-3232"><span class="mn" id="MathJax-Span-3233" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3234"><span class="msubsup" id="MathJax-Span-3235"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3236" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3237"><span class="mrow" id="MathJax-Span-3238"><span class="mn" id="MathJax-Span-3239" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3240"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3241" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3242"><span class="mrow" id="MathJax-Span-3243"><span class="mn" id="MathJax-Span-3244" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-3245" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3246" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3247" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3248"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3249" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3250"><span class="mrow" id="MathJax-Span-3251"><span class="mn" id="MathJax-Span-3252" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3253"><span class="mrow" id="MathJax-Span-3254"><span class="mn" id="MathJax-Span-3255" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3256"><span class="msubsup" id="MathJax-Span-3257"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3258" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3259"><span class="mrow" id="MathJax-Span-3260"><span class="mn" id="MathJax-Span-3261" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3262"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3263" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3264"><span class="mrow" id="MathJax-Span-3265"><span class="mn" id="MathJax-Span-3266" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-3267" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3268" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3269" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3270"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3271" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3272"><span class="mrow" id="MathJax-Span-3273"><span class="mn" id="MathJax-Span-3274" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3275"><span class="mrow" id="MathJax-Span-3276"><span class="mn" id="MathJax-Span-3277" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3278"><span class="msubsup" id="MathJax-Span-3279"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3280" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3281"><span class="mrow" id="MathJax-Span-3282"><span class="mn" id="MathJax-Span-3283" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3284"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3285" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3286"><span class="mrow" id="MathJax-Span-3287"><span class="mn" id="MathJax-Span-3288" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-3289" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3290" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3291" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3292"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3293" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3294"><span class="mrow" id="MathJax-Span-3295"><span class="mn" id="MathJax-Span-3296" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3297"><span class="mrow" id="MathJax-Span-3298"><span class="mn" id="MathJax-Span-3299" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3300"><span class="msubsup" id="MathJax-Span-3301"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3302" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3303"><span class="mrow" id="MathJax-Span-3304"><span class="mn" id="MathJax-Span-3305" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3306"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3307" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3308"><span class="mrow" id="MathJax-Span-3309"><span class="mn" id="MathJax-Span-3310" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-3311" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3312" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3313" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.015em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3314"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3315" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3316"><span class="mrow" id="MathJax-Span-3317"><span class="mn" id="MathJax-Span-3318" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3319"><span class="mrow" id="MathJax-Span-3320"><span class="mn" id="MathJax-Span-3321" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3322"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3323" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3324"><span class="mrow" id="MathJax-Span-3325"><span class="mn" id="MathJax-Span-3326" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.01em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.015em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3327" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3328" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(2.979em, 1000.9em, 4.467em, -999.997em); top: -4.818em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3329"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3330" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -4.342em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3331"><span class="mrow" id="MathJax-Span-3332"><span class="mn" id="MathJax-Span-3333" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.336em, 1000.42em, 4.17em, -999.997em); top: -3.688em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3334"><span class="mrow" id="MathJax-Span-3335"><span class="mn" id="MathJax-Span-3336" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3337"><span class="msubsup" id="MathJax-Span-3338"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3339" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3340"><span class="mrow" id="MathJax-Span-3341"><span class="mn" id="MathJax-Span-3342" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3343"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3344" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3345"><span class="mrow" id="MathJax-Span-3346"><span class="mn" id="MathJax-Span-3347" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-3348" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.205em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -1.068em; border-left: 0px solid; width: 0px; height: 3.146em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mo>−</mo><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>−</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msubsup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msubsup><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-5">- s_{3}^{2} + \frac{s_{3}}{\theta_{2 0}} - \frac{s_{3}^{2}}{s_{2} \theta_{2 0}} + \frac{s_{3}^{2}}{s_{2} \theta_{1 0}} - \frac{s_{3}^{2}}{s_{1} \theta_{1 0}} + \frac{s_{3}^{2}}{s_{1} \theta_{0 0}} + \frac{s_{3}^{2}}{s_{0}} - \frac{s_{3}^{2}}{s_{0} \theta_{0 0}}</script>
</div>

</div>

</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[10]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">beta</span> <span class="o">=</span> <span class="n">symbols</span><span class="p">(</span><span class="s1">'beta'</span><span class="p">)</span>
<span class="n">p_0</span> <span class="o">=</span> <span class="n">symbols</span><span class="p">(</span><span class="s1">'p_0'</span><span class="p">)</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'output = p_0 * beta * expand(s_simp / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="n">depth_of_tree</span><span class="p">)</span>
<span class="n">output</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt output_prompt">Out[10]:</div>




<div class="output_latex output_subarea output_execute_result">
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-6-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mi&gt;&amp;#x03B2;&lt;/mi&gt;&lt;msub&gt;&lt;mi&gt;p&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;2&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-3349" style="width: 34.943em; display: inline-block;"><span style="display: inline-block; position: relative; width: 29.11em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(1.789em, 1028.87em, 4.527em, -999.997em); top: -3.39em; left: 0em;"><span class="mrow" id="MathJax-Span-3350"><span class="mi" id="MathJax-Span-3351" style="font-family: MathJax_Math-italic;">β<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.003em;"></span></span><span class="msubsup" id="MathJax-Span-3352"><span style="display: inline-block; position: relative; width: 0.955em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.48em, 4.348em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3353" style="font-family: MathJax_Math-italic;">p</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3354"><span class="mrow" id="MathJax-Span-3355"><span class="mn" id="MathJax-Span-3356" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-3357" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-3358" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-3359"><span class="mo" id="MathJax-Span-3360" style="font-family: MathJax_Main;">−</span><span class="msubsup" id="MathJax-Span-3361"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3362" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3363"><span class="mrow" id="MathJax-Span-3364"><span class="mn" id="MathJax-Span-3365" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-3366" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3367" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.372em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.158em, 1000.42em, 4.17em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.235em;"><span class="mn" id="MathJax-Span-3368" style="font-family: MathJax_Main;">1</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1001.25em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.652em;"><span class="msubsup" id="MathJax-Span-3369"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3370" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3371"><span class="mrow" id="MathJax-Span-3372"><span class="mn" id="MathJax-Span-3373" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-3374" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.37em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.372em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3375" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3376" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3377"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3378" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3379"><span class="mrow" id="MathJax-Span-3380"><span class="mn" id="MathJax-Span-3381" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3382"><span class="msubsup" id="MathJax-Span-3383"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3384" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3385"><span class="mrow" id="MathJax-Span-3386"><span class="mn" id="MathJax-Span-3387" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3388"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3389" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3390"><span class="mrow" id="MathJax-Span-3391"><span class="mn" id="MathJax-Span-3392" style="font-size: 70.7%; font-family: MathJax_Main;">2</span><span class="mn" id="MathJax-Span-3393" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3394" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3395" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3396"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3397" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3398"><span class="mrow" id="MathJax-Span-3399"><span class="mn" id="MathJax-Span-3400" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3401"><span class="msubsup" id="MathJax-Span-3402"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3403" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3404"><span class="mrow" id="MathJax-Span-3405"><span class="mn" id="MathJax-Span-3406" style="font-size: 70.7%; font-family: MathJax_Main;">2</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3407"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3408" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3409"><span class="mrow" id="MathJax-Span-3410"><span class="mn" id="MathJax-Span-3411" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-3412" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3413" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3414" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3415"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3416" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3417"><span class="mrow" id="MathJax-Span-3418"><span class="mn" id="MathJax-Span-3419" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3420"><span class="msubsup" id="MathJax-Span-3421"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3422" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3423"><span class="mrow" id="MathJax-Span-3424"><span class="mn" id="MathJax-Span-3425" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3426"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3427" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3428"><span class="mrow" id="MathJax-Span-3429"><span class="mn" id="MathJax-Span-3430" style="font-size: 70.7%; font-family: MathJax_Main;">1</span><span class="mn" id="MathJax-Span-3431" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3432" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3433" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3434"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3435" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3436"><span class="mrow" id="MathJax-Span-3437"><span class="mn" id="MathJax-Span-3438" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3439"><span class="msubsup" id="MathJax-Span-3440"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3441" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3442"><span class="mrow" id="MathJax-Span-3443"><span class="mn" id="MathJax-Span-3444" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3445"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3446" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3447"><span class="mrow" id="MathJax-Span-3448"><span class="mn" id="MathJax-Span-3449" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-3450" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3451" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3452" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.015em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3453"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3454" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3455"><span class="mrow" id="MathJax-Span-3456"><span class="mn" id="MathJax-Span-3457" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3458"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3459" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3460"><span class="mrow" id="MathJax-Span-3461"><span class="mn" id="MathJax-Span-3462" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.01em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.015em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3463" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3464" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3465"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3466" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3467"><span class="mrow" id="MathJax-Span-3468"><span class="mn" id="MathJax-Span-3469" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3470"><span class="msubsup" id="MathJax-Span-3471"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3472" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3473"><span class="mrow" id="MathJax-Span-3474"><span class="mn" id="MathJax-Span-3475" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3476"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3477" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3478"><span class="mrow" id="MathJax-Span-3479"><span class="mn" id="MathJax-Span-3480" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-3481" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-3482" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span></span><span style="display: inline-block; width: 0px; height: 3.396em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -1.211em; border-left: 0px solid; width: 0px; height: 3.004em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>β</mi><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><mrow><mo>(</mo><mrow><mo>−</mo><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mo>+</mo><mfrac><mn>1</mn><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>2</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac></mrow><mo>)</mo></mrow></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-6">\beta p_{0} \left(- s_{3} + \frac{1}{\theta_{2 0}} - \frac{s_{3}}{s_{2} \theta_{2 0}} + \frac{s_{3}}{s_{2} \theta_{1 0}} - \frac{s_{3}}{s_{1} \theta_{1 0}} + \frac{s_{3}}{s_{1} \theta_{0 0}} + \frac{s_{3}}{s_{0}} - \frac{s_{3}}{s_{0} \theta_{0 0}}\right)</script>
</div>

</div>

</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[73]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1">#### Deriving the cross-elasticity formular for vehicles in the same make_region (second to the lowest level)</span>
<span class="n">level_of_similarity</span> <span class="o">=</span> <span class="mi">0</span> <span class="c1"># starts from 0</span>
<span class="n">diff_node</span> <span class="o">=</span> <span class="n">total_end_nodes</span> <span class="o">//</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="n">level_of_similarity</span> <span class="o">-</span> <span class="mi">1</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'s_diff_v_x = diff(s, v_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="n">diff_node</span><span class="p">)</span>
<span class="n">s_v_0</span> <span class="o">=</span> <span class="n">s</span>
<span class="n">s_simp</span> <span class="o">=</span> <span class="n">s_diff_v_x</span><span class="o">.</span><span class="n">subs</span><span class="p">(</span><span class="n">sum_exp_0_0</span> <span class="o">/</span> <span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">sum_exp_0_0</span><span class="p">),</span> <span class="n">s_0</span><span class="p">)</span>
<span class="n">s_v_0</span> <span class="o">=</span> <span class="n">s_v_0</span><span class="o">.</span><span class="n">subs</span><span class="p">(</span><span class="n">sum_exp_0_0</span> <span class="o">/</span> <span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">sum_exp_0_0</span><span class="p">),</span> <span class="n">s_0</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">level_of_similarity</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_simp = s_simp.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_0 / sum_</span><span class="si">%d</span><span class="s1">_0, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_v_0 = s_v_0.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_0 / sum_</span><span class="si">%d</span><span class="s1">_0, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">level_of_similarity</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_simp = s_simp.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> / sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span> <span class="o">//</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="n">i</span><span class="o">-</span><span class="mi">2</span><span class="p">),</span>
                                                                           <span class="n">i</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span> <span class="o">//</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="n">i</span><span class="o">-</span><span class="mi">1</span><span class="p">),</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_v_0 = s_v_0.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1"> / sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span> <span class="o">//</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="n">i</span><span class="o">-</span><span class="mi">2</span><span class="p">),</span>
                                                                           <span class="n">i</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span> <span class="o">//</span> <span class="n">number_splits_non_end_node</span> <span class="o">**</span> <span class="p">(</span><span class="n">depth_of_tree</span><span class="o">-</span><span class="n">i</span><span class="o">-</span><span class="mi">1</span><span class="p">),</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'s_simp = s_simp.subs(exp(v_</span><span class="si">%d</span><span class="s1">/theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">) / sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">diff_node</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">,</span>
                                                                           <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="p">,</span>
                                                                              <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">))</span>
<span class="n">exec</span><span class="p">(</span><span class="s1">'s_v_0 = s_v_0.subs(exp(v_</span><span class="si">%d</span><span class="s1">/theta_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">) / sum_</span><span class="si">%d</span><span class="s1">_</span><span class="si">%d</span><span class="s1">, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">diff_node</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">,</span>
                                                                           <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">,</span> <span class="n">diff_node</span> <span class="o">//</span> <span class="n">number_end_group</span><span class="p">,</span> <span class="n">depth_of_tree</span><span class="p">,</span>
                                                                              <span class="n">depth_of_tree</span><span class="o">-</span><span class="mi">1</span><span class="p">))</span>
<span class="n">beta</span> <span class="o">=</span> <span class="n">Symbol</span><span class="p">(</span><span class="s1">'beta'</span><span class="p">)</span>
<span class="n">p_1</span> <span class="o">=</span> <span class="n">Symbol</span><span class="p">(</span><span class="s1">'p_1'</span><span class="p">)</span>
<span class="n">beta</span> <span class="o">*</span> <span class="n">p_1</span> <span class="o">*</span> <span class="n">simplify</span><span class="p">(</span><span class="n">s_simp</span> <span class="o">/</span> <span class="n">s_v_0</span><span class="p">)</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt output_prompt">Out[73]:</div>




<div class="output_latex output_subarea output_execute_result">
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-7-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mi&gt;&amp;#x03B2;&lt;/mi&gt;&lt;msub&gt;&lt;mi&gt;p&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;mo&gt;(&lt;/mo&gt;&lt;mrow&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mo&gt;+&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mfrac&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mfrac&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;mrow&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;&amp;#x03B8;&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;mn&gt;0&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/mrow&gt;&lt;/mfrac&gt;&lt;/mrow&gt;&lt;mo&gt;)&lt;/mo&gt;&lt;/mrow&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-3483" style="width: 13.396em; display: inline-block;"><span style="display: inline-block; position: relative; width: 11.134em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(1.789em, 1010.9em, 4.527em, -999.997em); top: -3.39em; left: 0em;"><span class="mrow" id="MathJax-Span-3484"><span class="mi" id="MathJax-Span-3485" style="font-family: MathJax_Math-italic;">β<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.003em;"></span></span><span class="msubsup" id="MathJax-Span-3486"><span style="display: inline-block; position: relative; width: 0.955em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.48em, 4.348em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3487" style="font-family: MathJax_Math-italic;">p</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3488"><span class="mrow" id="MathJax-Span-3489"><span class="mn" id="MathJax-Span-3490" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mrow" id="MathJax-Span-3491" style="padding-left: 0.182em;"><span class="mo" id="MathJax-Span-3492" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">(</span></span><span class="mrow" id="MathJax-Span-3493"><span class="mo" id="MathJax-Span-3494" style="font-family: MathJax_Main;">−</span><span class="msubsup" id="MathJax-Span-3495"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3496" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3497"><span class="mrow" id="MathJax-Span-3498"><span class="mn" id="MathJax-Span-3499" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="mo" id="MathJax-Span-3500" style="font-family: MathJax_Main; padding-left: 0.241em;">+</span><span class="mfrac" id="MathJax-Span-3501" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 1.015em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3502"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3503" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3504"><span class="mrow" id="MathJax-Span-3505"><span class="mn" id="MathJax-Span-3506" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3507"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3508" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3509"><span class="mrow" id="MathJax-Span-3510"><span class="mn" id="MathJax-Span-3511" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1001.01em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 1.015em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span><span class="mo" id="MathJax-Span-3512" style="font-family: MathJax_Main; padding-left: 0.241em;">−</span><span class="mfrac" id="MathJax-Span-3513" style="padding-left: 0.241em;"><span style="display: inline-block; position: relative; width: 2.265em; height: 0px; margin-right: 0.122em; margin-left: 0.122em;"><span style="position: absolute; clip: rect(3.396em, 1000.9em, 4.348em, -999.997em); top: -4.64em; left: 50%; margin-left: -0.473em;"><span class="msubsup" id="MathJax-Span-3514"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3515" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3516"><span class="mrow" id="MathJax-Span-3517"><span class="mn" id="MathJax-Span-3518" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(3.098em, 1002.15em, 4.348em, -999.997em); top: -3.271em; left: 50%; margin-left: -1.068em;"><span class="mrow" id="MathJax-Span-3519"><span class="msubsup" id="MathJax-Span-3520"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3521" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3522"><span class="mrow" id="MathJax-Span-3523"><span class="mn" id="MathJax-Span-3524" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3525"><span style="display: inline-block; position: relative; width: 1.253em; height: 0px;"><span style="position: absolute; clip: rect(3.098em, 1000.48em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3526" style="font-family: MathJax_Math-italic;">θ</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3527"><span class="mrow" id="MathJax-Span-3528"><span class="mn" id="MathJax-Span-3529" style="font-size: 70.7%; font-family: MathJax_Main;">0</span><span class="mn" id="MathJax-Span-3530" style="font-size: 70.7%; font-family: MathJax_Main;">0</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; clip: rect(0.836em, 1002.26em, 1.253em, -999.997em); top: -1.307em; left: 0em;"><span style="display: inline-block; overflow: hidden; vertical-align: 0em; border-top: 1.3px solid; width: 2.265em; height: 0px;"></span><span style="display: inline-block; width: 0px; height: 1.074em;"></span></span></span></span></span><span class="mo" id="MathJax-Span-3531" style="vertical-align: 0em;"><span style="font-family: MathJax_Size3;">)</span></span></span></span><span style="display: inline-block; width: 0px; height: 3.396em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -1.211em; border-left: 0px solid; width: 0px; height: 3.004em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>β</mi><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><mrow><mo>(</mo><mrow><mo>−</mo><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mo>+</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub></mfrac><mo>−</mo><mfrac><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub><mrow><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn></mrow></msub><msub><mi>θ</mi><mrow class="MJX-TeXAtom-ORD"><mn>0</mn><mn>0</mn></mrow></msub></mrow></mfrac></mrow><mo>)</mo></mrow></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-7">\beta p_{1} \left(- s_{3} + \frac{s_{3}}{s_{0}} - \frac{s_{3}}{s_{0} \theta_{0 0}}\right)</script>
</div>

</div>

</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[79]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="c1"># Cross-elasticity with outside option</span>
<span class="n">outside_s</span> <span class="o">=</span> <span class="mi">1</span> <span class="o">-</span> <span class="n">sum_exp_0_0</span> <span class="o">/</span> <span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">sum_exp_0_0</span><span class="p">)</span>
<span class="n">s_diff_v_0</span> <span class="o">=</span> <span class="n">diff</span><span class="p">(</span><span class="n">outside_s</span><span class="p">,</span> <span class="n">v_0</span><span class="p">)</span>
<span class="n">s_simp</span> <span class="o">=</span> <span class="n">s_diff_v_0</span> <span class="o">/</span> <span class="n">outside_s</span>
<span class="n">s_simp</span> <span class="o">=</span> <span class="n">s_simp</span><span class="o">.</span><span class="n">subs</span><span class="p">(</span><span class="n">sum_exp_0_0</span> <span class="o">/</span> <span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">sum_exp_0_0</span><span class="p">),</span> <span class="n">s_0</span><span class="p">)</span>
<span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">depth_of_tree</span><span class="p">):</span>
    <span class="n">exec</span><span class="p">(</span><span class="s1">'s_simp = s_simp.subs(sum_exp_</span><span class="si">%d</span><span class="s1">_0 / sum_</span><span class="si">%d</span><span class="s1">_0, s_</span><span class="si">%d</span><span class="s1"> / s_</span><span class="si">%d</span><span class="s1">)'</span> <span class="o">%</span> <span class="p">(</span><span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">i</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span> <span class="n">i</span><span class="p">))</span>
<span class="n">beta</span> <span class="o">=</span> <span class="n">Symbol</span><span class="p">(</span><span class="s1">'beta'</span><span class="p">)</span>
<span class="n">p_1</span> <span class="o">=</span> <span class="n">Symbol</span><span class="p">(</span><span class="s1">'p_1'</span><span class="p">)</span>
<span class="n">beta</span> <span class="o">*</span> <span class="n">p_1</span> <span class="o">*</span> <span class="n">simplify</span><span class="p">(</span><span class="n">s_simp</span><span class="p">)</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area">

<div class="prompt output_prompt">Out[79]:</div>




<div class="output_latex output_subarea output_execute_result">
<span class="MathJax_Preview" style="color: inherit; display: none;"></span><div class="MathJax_Display" style="text-align: center;"><span class="MathJax" id="MathJax-Element-8-Frame" tabindex="0" data-mathml="&lt;math xmlns=&quot;http://www.w3.org/1998/Math/MathML&quot; display=&quot;block&quot;&gt;&lt;mo&gt;&amp;#x2212;&lt;/mo&gt;&lt;mi&gt;&amp;#x03B2;&lt;/mi&gt;&lt;msub&gt;&lt;mi&gt;p&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;1&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;msub&gt;&lt;mi&gt;s&lt;/mi&gt;&lt;mrow class=&quot;MJX-TeXAtom-ORD&quot;&gt;&lt;mn&gt;3&lt;/mn&gt;&lt;/mrow&gt;&lt;/msub&gt;&lt;/math&gt;" role="presentation" style="text-align: center; position: relative;"><nobr aria-hidden="true"><span class="math" id="MathJax-Span-3532" style="width: 3.872em; display: inline-block;"><span style="display: inline-block; position: relative; width: 3.217em; height: 0px; font-size: 120%;"><span style="position: absolute; clip: rect(1.313em, 1003.22em, 2.562em, -999.997em); top: -2.199em; left: 0em;"><span class="mrow" id="MathJax-Span-3533"><span class="mo" id="MathJax-Span-3534" style="font-family: MathJax_Main;">−</span><span class="mi" id="MathJax-Span-3535" style="font-family: MathJax_Math-italic;">β<span style="display: inline-block; overflow: hidden; height: 1px; width: 0.003em;"></span></span><span class="msubsup" id="MathJax-Span-3536"><span style="display: inline-block; position: relative; width: 0.955em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.48em, 4.348em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3537" style="font-family: MathJax_Math-italic;">p</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3538"><span class="mrow" id="MathJax-Span-3539"><span class="mn" id="MathJax-Span-3540" style="font-size: 70.7%; font-family: MathJax_Main;">1</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span><span class="msubsup" id="MathJax-Span-3541"><span style="display: inline-block; position: relative; width: 0.896em; height: 0px;"><span style="position: absolute; clip: rect(3.396em, 1000.42em, 4.17em, -999.997em); top: -3.985em; left: 0em;"><span class="mi" id="MathJax-Span-3542" style="font-family: MathJax_Math-italic;">s</span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span><span style="position: absolute; top: -3.807em; left: 0.479em;"><span class="texatom" id="MathJax-Span-3543"><span class="mrow" id="MathJax-Span-3544"><span class="mn" id="MathJax-Span-3545" style="font-size: 70.7%; font-family: MathJax_Main;">3</span></span></span><span style="display: inline-block; width: 0px; height: 3.991em;"></span></span></span></span></span><span style="display: inline-block; width: 0px; height: 2.205em;"></span></span></span><span style="display: inline-block; overflow: hidden; vertical-align: -0.282em; border-left: 0px solid; width: 0px; height: 1.218em;"></span></span></nobr><span class="MJX_Assistive_MathML MJX_Assistive_MathML_Block" role="presentation"><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mo>−</mo><mi>β</mi><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mn>1</mn></mrow></msub><msub><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mn>3</mn></mrow></msub></math></span></span></div><script type="math/tex; mode=display" id="MathJax-Element-8">- \beta p_{1} s_{3}</script>
</div>

</div>

</div>
</div>

</div>
    </div>
  </div>


 



<div style="position: absolute; width: 0px; height: 0px; overflow: hidden; padding: 0px; border: 0px; margin: 0px;"><div id="MathJax_Font_Test" style="position: absolute; visibility: hidden; top: 0px; left: 0px; width: auto; padding: 0px; border: 0px; margin: 0px; white-space: nowrap; text-align: left; text-indent: 0px; text-transform: none; line-height: normal; letter-spacing: normal; word-spacing: normal; font-size: 40px; font-weight: normal; font-style: normal; font-family: MathJax_Size4, sans-serif;"></div></div></body></html>