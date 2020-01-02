//copy to Clipboard
function fallbackCopyTextToClipboard(text) {
  var textArea = document.createElement("textarea");
  textArea.value = text;
  document.body.appendChild(textArea);
  textArea.focus();
  textArea.select();

  try {
    var successful = document.execCommand('copy');
    var msg = successful ? 'successful' : 'unsuccessful';
    console.log('Fallback: Copying text command was ' + msg);
  } catch (err) {
    console.error('Fallback: Oops, unable to copy', err);
  }

  document.body.removeChild(textArea);
}
function copyTextToClipboard(text) {
  if (!navigator.clipboard) {
    fallbackCopyTextToClipboard(text);
    return;
  }
  navigator.clipboard.writeText(text).then(function() {
    console.log('Async: Copying to clipboard was successful!');
  }, function(err) {
    console.error('Async: Could not copy text: ', err);
  });
}

var data = source.data;
var copytext = 'x,y,e\n';
for (var i = 0; i < data['n'].length; i++) {
    var currRow = [data['n'][i].toString(),
                   data['yt'][i].toString(),
                   data['et'][i].toString().concat('\n')];
    var joined = currRow.join();
    copytext = copytext.concat(joined);
}

copyTextToClipboard(copytext)
